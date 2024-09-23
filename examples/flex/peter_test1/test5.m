% Parameters

n_objs = 2; % number of objects
source_loc = [-3;0]; % location of test source
source_loc = [2.5;2.5];
zk = 2*pi;  % our k (wave number)
nu = 1/3; % Poisson ratio
cparams = [];
cparams.maxchunklen = 0.5 / zk; % max chunk size
narms = 5; % number of arms on the starfish

s = 5; % spacing
box_size = s*ceil(sqrt(n_objs))-s;
[xs, ys] = meshgrid(0:s:box_size,0:s:box_size);
centers = [xs(:) ys(:)];

temp = [centers, abs(centers(:,1) - centers(:,2))];
temp = sortrows(temp, 3);
centers = temp(1:n_objs,1:2);

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
for i = 1:n_objs
    chnkr = chunkerfunc(@(t) starfish(t,narms), cparams);
    chnkr = chnkr.move([0;0],centers(i,:)',2*i,1);
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

for ichk = 1:n_objs

    indst = iinds(ichk) + (1:ninds(ichk)) - 1;
    targs = [];
    targs.r = sout{ichk}.sk_in_r;
    targs.d = sout{ichk}.sk_in_d;
    targs.n = sout{ichk}.sk_in_n;
    targs.kappa = sout{ichk}.kappa;
    [data] = flex_pt(zk,coefs,source_loc,targs);
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

xs = -s/2:s/100:box_size+s/2;                                    % generate some targets
ys = -s/2:s/100:box_size+s/2; 
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

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

utarg = zeros(na, 1);
utarg(out) = utot-trueval;
utarg(~out) = NaN;

utarg = reshape(utarg,sz);
h = pcolor(X,Y,real(utarg));
set(h,'EdgeColor','None'); hold on;
title("BIE solution")
for ii=1:n_objs
    plot(chnkrs{ii},'w-','LineWidth',2); hold on;
end
colorbar


%abs(utot-( -0.003566175236184 + 0.001269459882790i))/abs( -0.003566175236184 + 0.001269459882790i)


