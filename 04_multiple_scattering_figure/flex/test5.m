% Parameters

n_objs = 101; % number of objects
source_loc = [-3;0]; % location of test source
source_loc = [2.5;2.5];
zk = 6;  % our k (wave number)
nu = 1/3; % Poisson ratio
cparams = [];
cparams.maxchunklen = 4 / zk; % max chunk size

s = 5; % spacing
box_size = s*ceil(sqrt(n_objs))-s;
[xs, ys] = meshgrid(0:s:box_size,0:s:box_size);
centers = [xs(:) ys(:)];
centers = centers + 0.*(rand(size(centers))-0.5);
temp = [centers, abs(centers(:,1) - centers(:,2))];
temp = sortrows(temp, 3);
centers = temp(1:n_objs,1:2);
plot(centers(:,1),centers(:,2),'o');

% Generating chnkr objects

npxy = 200;
pxyr = 2;
pxypts1 = pxyr*[cos(2*pi*(0:(npxy-1))/npxy);sin(2*pi*(0:(npxy-1))/npxy)];
pxypts2 = (1.2)*pxyr*[cos(2*pi*(0:(npxy-1))/npxy);sin(2*pi*(0:(npxy-1))/npxy)];
pxypts = [pxypts1,pxypts2];

coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate first part', coefs);                          % build the desired kernel
fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
fkern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);                    % first part K21 coupled with hilbert transforms

hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
double = @(s,t) chnk.lap2d.kern(s,t, 'd');

fkern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part', coefs);             % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)
fkern4 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 second part', coefs);            % kernels in K21 needs to multiply by curvature
fkern5 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);           % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature
fkern6 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K22 second part', coefs);            % kernels in K22 needs to multiply by curvature
fkern1_trans = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert unsubtract', coefs);  % K11 without the hilbert subtraction
fkern3_trans = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part unsubtract', coefs);             % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

kerns = {};
kerns{1} = fkern;
kerns{2} = fkern1;
kerns{3} = fkern2;
kerns{4} = hilbert;
kerns{5} = double;
kerns{6} = fkern3;
kerns{7} = fkern4;
kerns{8} = fkern5;
kerns{9} = fkern6;
kerns{10}= fkern1_trans;
kerns{11}= fkern3_trans;

kerns{12}= ikern1;
kerns{13}= ikern2;
kerns{14}= ikern3;

tgen = tic;

sout = {};
%chnkrs = createArray(n_objs,1,'chunker');
chnkrs = {};
narm = datasample([3 4 5 6], n_objs, 'Weights',[0.3 0.25 0.25 0.25]);

for i = 1:n_objs
    %narm = round(rand(1)*4 - 0.5)+3; % between 3 and 6 and change amp 
    amp = 0.15 + 0.3*rand;

    if (i == 1) || (i == 12) || (i == 13) 
        rot = -pi/12;
        narm(i) = 3;
    else 
        rot = 2*pi*rand;
    end

    chnkr = chunkerfunc(@(t) starfish(t,narm(i),amp,[0,0],0, 1), cparams);
    chnkr = chnkr.move([0;0],centers(i,:)',rot,1);
    plot(chnkr,'k-','LineWidth',2); hold on;
    drawnow
    %chnkrs(i) = chnkr;
    sout{i} = skel_chnkr(chnkr,nu,zk,coefs,kerns,pxypts);
    chnkrs{i} = chnkr;
end

t2 = toc(tgen);
fprintf('%5.2f s : time for skeletonizing \n',t2)

iinds = zeros(n_objs,1);
ninds = zeros(n_objs,1);

icount = 1;
for ii=1:n_objs
    iinds(ii) = icount;
    ninds(ii) = 3*size(sout{ii}.sk_ou_r,2);
    icount = icount + ninds(ii);
end
ntot = icount - 1;

%%
%%%%%%%%%%%% generate the system matrix

tsyscrea = tic;

sys = eye(ntot,ntot);

for itar=1:n_objs

    tar = [];
    tar.r = sout{itar}.sk_in_r;
    tar.n = sout{itar}.sk_in_n;
    tar.d = sout{itar}.sk_in_d;
    tar.kappa = sout{itar}.kappa;

    indst = iinds(itar) + (1:ninds(itar)) - 1;

    for isrc = 1:n_objs
        
        if (isrc ~= itar)

        src = [];
        src.r = sout{isrc}.sk_ou_r;
        src.n = sout{isrc}.sk_ou_n;
        src.d = sout{isrc}.sk_ou_d;
        src.wts = sout{isrc}.sk_ou_w;

        indss = iinds(isrc) + (1:ninds(isrc)) - 1;

        tmat = get_trans_mat(src,tar,kerns);
        tmat = sout{itar}.scat_red*tmat;

        sys(indst,indss) = tmat;
        
        end

    end
end

t2 = toc(tsyscrea);

fprintf('%5.2f s : time for full matrix generation \n',t2)

%%
%%%%%%%%%%%%%%%%%%% solve

rhs = zeros(ntot,1);
theta = pi/4;
dir = [cos(theta);sin(theta)];
for ichk = 1:n_objs

    indst = iinds(ichk) + (1:ninds(ichk)) - 1;
    targs = [];
    targs.r = sout{ichk}.sk_in_r;
    targs.d = sout{ichk}.sk_in_d;
    targs.n = sout{ichk}.sk_in_n;
    targs.kappa = sout{ichk}.kappa;
    %[data] = flex_pt(zk,coefs,source_loc,targs);
    [data] = flex_wv(zk,coefs,dir,targs);
    tdens = sout{ichk}.scat_red*data;
    rhs(indst,1) = tdens;

end



%%%%%%%% now solve

tsolve = tic;
dens_red = sys\rhs;
t2 = toc(tsolve);
fprintf('%5.2f s : time for solve \n',t2)

%% post process


targ   = [];

xs = -s:s/80:box_size+s;                                    % generate some targets
ys = -s:s/80:box_size+s; 
[X,Y] = meshgrid(xs, ys);
sz = size(X);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

in = 0;
for i = 1:n_objs
    in = in + chunkerinterior(chnkrs{i}, targets); 
end
out = ~in; 

targ.r = targets(:,find(out));

utot = 0;

for ichk = 1:n_objs

    indst = iinds(ichk) + (1:ninds(ichk)) - 1;
    tdens = dens_red(indst);

    v1 = tdens(1:3:end);
    v2 = tdens(2:3:end);
    v3 = tdens(3:3:end);

    ch1src = [];
    ch1src.r = sout{ichk}.sk_ou_r;
    ch1src.n = sout{ichk}.sk_ou_n;
    ch1src.d = sout{ichk}.sk_ou_d;
    ch1src.w = sout{ichk}.sk_ou_w;
    ik1b = ikern1(ch1src,targ)*(v1.*ch1src.w);
    ik2b = ikern2(ch1src,targ)*(v3.*ch1src.w);
    ik3b = ikern3(ch1src,targ)*(v2.*ch1src.w);

    utot = utot + ik1b+ik2b+ik3b;

end

[val, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, source_loc, targ.r);
[valK, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), source_loc, targ.r);

%%
%trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;
trueval = exp(1i*zk*(targ.r(1,:)*dir(1)+targ.r(2,:)*dir(2)));
utarg = zeros(na, 1);
utarg(out) = -utot+trueval.';
utarg(~out) = NaN;
utarg = reshape(utarg,sz);

%%

figure(1)
for ii=1:n_objs
   plot(chnkrs{ii},'k-','LineWidth',2); hold on;
   %scatter(sout{ii}.pxypts(1,:),sout{ii}.pxypts(2,:),1,'color','red')
end
hold on
xlim([-4, max(X(:))])
ylim([-4, max(Y(:))])
caxis([0,3])
hold on 
plot([40 40+11],[0 11])
colorbar

figure(2)
h = pcolor(X,Y,abs(utarg));
%h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
axis off
xlim([-4, max(X(:))])
ylim([-4, max(Y(:))])
caxis([0,3])