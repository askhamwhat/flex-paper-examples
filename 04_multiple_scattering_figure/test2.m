% Parameters

n_objs = 1; % number of objects
source_loc = [-3;0]; % location of test source

zk = 1;  % our k (wave number)
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

%chnkrs = createArray(n_objs,1,'chunker');
chnkrs = [];
for i = 1:n_objs
    chnkr = chunkerfunc(@(t) starfish(t,narms), cparams);
    chnkr = chnkr.move([0;0],centers(i,:)',0,1);
    %chnkrs(i) = chnkr;
    chnkrs = [chnkrs,chnkr];
end

clf
figure(1)
plot(chnkrs, '-x'); hold on
quiver(chnkrs); hold on
axis equal
drawnow


% Defining kernels

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


chnkr = chnkrs;

% Building lefthand side

start = tic;


src_chnkr = chnkr;
opts2 = [];
opts2.sing = 'pv';
H = chunkermat(src_chnkr, hilbert, opts2);
wts = src_chnkr.wts(:)';



sysmat = chunkermat(chnkr,fkern, opts);
sysmat1 = chunkermat(chnkr, fkern1, opts);
sysmat2 = chunkermat(chnkr, fkern2, opts);

K21 = chunkermat(chnkr, fkern3, opts);
K21second = chunkermat(chnkr, fkern4, opts);
K21hilbert = chunkermat(chnkr, fkern5, opts);
K22second = chunkermat(chnkr, fkern6, opts);
D = chunkermat(chnkr, double, opts);                                            % Assemble hilbert transforms

kappa = signed_curvature(chnkr);
kappa = kappa(:)';

hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa.').*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + (kappa.').*(K22second);

A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
M = kron(eye(chnkr.npt), A);
Sinv =  M + sysmat ;

lhs = Sinv;

disp(['S_',num2str(i)])





    t1 = toc(start);
    fprintf('%5.2f s : time to assemble matrix\n',t1)


    % Evaluating righthand side(s)

    start = tic;
    rhs = [];

    [~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, source_loc, chnkr.r);
    [~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), source_loc, chnkr.r);

    nx = chnkr.n(1,:).';
    ny = chnkr.n(2,:).';

    dx = chnkr.d(1,:).';
    dy = chnkr.d(2,:).';

    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds);                                                                       % normalization
    tauy = (dy./ds);

    kappa = signed_curvature(chnkr);
    kappa = kappa(:)';

    firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))-...
        1/(2*zk^2).*(hessK(:, :, 1).*(nx.*nx) + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*(ny.*ny))+...
        coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))-...
        coefs(1)/(2*zk^2).*(hessK(:, :, 1).*(taux.*taux) + hessK(:, :, 2).*(2*taux.*tauy) + ...
        hessK(:, :, 3).*(tauy.*tauy));

    secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
        third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) - ...
        1./(2*zk^2).*(thirdK(:, :, 1).*(nx.*nx.*nx) + thirdK(:, :, 2).*(3*nx.*nx.*ny)+...
        thirdK(:, :, 3).*(3*nx.*ny.*ny) + thirdK(:, :, 4).*(ny.*ny.*ny)) +...
        (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:, :, 4).*(tauy.*tauy.*ny)) - ...
        (2-coefs(1))/(2*zk^2).*(thirdK(:, :, 1).*(taux.*taux.*nx) + thirdK(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        thirdK(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + thirdK(:, :, 4).*(tauy.*tauy.*ny)) +...
        (1-coefs(1))*diag(kappa)*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
        1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
        hessK(:, :, 3).*tauy.*tauy)-...
        (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
        1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));

    data = zeros(2*chnkr.npt, 1);
    data(1:2:2*chnkr.npt) = firstbc ;
    data(2:2:2*chnkr.npt) = secondbc;

    rhs = [rhs; data];


    t1 = toc(start);
    fprintf('%5.2f s : time to evaluate BCs \n',t1)

    % Solving linear system

    start = tic;

    sol = lhs\rhs;

    t1 = toc(start);
    fprintf('%5.2f s : time to solve linear system \n',t1)

    start1 = tic;



    % Plotting the solution

    xs = -s/2:s/20:box_size+s/2;                                    % generate some targets
    ys = -s/2:s/20:box_size+s/2;
    [X,Y] = meshgrid(xs, ys);
    targets = [X(:).'; Y(:).'];
    [~,na] = size(targets);

    start = tic;
    in = 0;
    in = in + chunkerinterior(chnkr, targets);
    out = ~in;
    t1 = toc(start);
    fprintf('%5.2f s : time to find interior points \n',t1)

    ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation
    ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
    ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

    %H = chunkermat(chnkr1, hilbert, opts2);                                              % Assemble hilbert transforms

    uscat = 0;
    count = 0;
    start1 = tic;

    rho1 = sol(count+1:2:count+2*chnkr.npt);
    rho2 = sol(count+2:2:count+2*chnkr.npt);
    count = count + 2*chnkrs(i).npt;

    H = chunkermat(chnkr, hilbert, opts2);                                              % Assemble hilbert transforms

    uscat =uscat + chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) +...
        chunkerkerneval(chnkr, ikern3,H*rho1, targets(:, out)) + ...
        chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));

ifplot = false;
if (ifplot)
    [val, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, source_loc, targets(:,out));
    [valK, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), source_loc, targets(:,out));

    trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

    true_sol = zeros(na, 1);
    utarg = zeros(na, 1);

    utarg(out) = uscat;
    true_sol(out) = trueval;

    utarg = reshape(utarg,size(X));
    true_sol = reshape(true_sol,size(X));

    t2 = toc(start1);
    fprintf('%5.2f s : time for kernel eval (for plotting)\n',t2)

    figure(2)
    tiledlayout(1,3)
    nexttile
    h = pcolor(X,Y,real(true_sol));
    set(h,'EdgeColor','None'); hold on;
    title("True solution")
    plot(chnkr,'w-','LineWidth',2);
    colorbar

    nexttile
    h = pcolor(X,Y,real(utarg));
    set(h,'EdgeColor','None'); hold on;
    title("BIE solution")
    plot(chnkr,'w-','LineWidth',2);
    colorbar

    nexttile
    utot = -utarg + true_sol;
    utot = reshape(utot,size(X));
    h = pcolor(X,Y,real(utot));
    set(h,'EdgeColor','None'); hold on;
    title("Log_{10} relative error (coupled)")
    plot(chnkr,'w-','LineWidth',2);
    colorbar

end

%%

npxy = 200;
pxyr = 2;
pxypts1 = pxyr*[cos(2*pi*(0:(npxy-1))/npxy);sin(2*pi*(0:(npxy-1))/npxy)];
pxypts2 = (1.2)*pxyr*[cos(2*pi*(0:(npxy-1))/npxy);sin(2*pi*(0:(npxy-1))/npxy)];
pxypts = [pxypts1,pxypts2];
T1 = chunkerkernevalmat(chnkr, ikern1, pxypts);
T2 = chunkerkernevalmat(chnkr, ikern3, pxypts);
T3 = chunkerkernevalmat(chnkr, ikern2, pxypts);
% 
% sz = size(T1);
% sz(2) = sz(2) + size(T2,2) + size(T3,2);
% T = zeros(sz);
% T(:,1:3:end) = T1;
T1 = T1;
% T(:,2:3:end) = T2;
T2 = T2;
% T(:,3:3:end) = T3;
T3 = T3;
% 
% tol = 1E-10;
% [sk,rd,skmat] = id(T,tol);
% tmat = zeros(size(sk,2),sz(2));
% tmat(:,rd) = skmat;
% tmat(:,sk) = eye(numel(sk));
T = [T1;T2;T3];
tol = 1E-10;
I = eye(size(T1,1));
[sk,~,~] = id(T,tol);
[sk1,~,~] = id(T1,tol);
[Q1,~] = qr(T1(:,sk1));
T2red = (I-Q1*Q1')*T2;
[sk2,~,~] = id(T2red,tol/norm(T2red,'fro'));
sk2tot = unique([sk1,sk2]);
[Q2,~] = qr([T1(:,sk2tot),T2(:,sk2tot)]);
T3red = (I-Q2*Q2')*T3;
[sk3,~,~] = id(T3red,tol/norm(T3red,'fro'));
sk = unique([sk1,sk2,sk3]);

sz = size(T1);
sz(2) = sz(2) + size(T2,2) + size(T3,2);
T = zeros(sz);
T(:,1:3:end) = T1;
T(:,2:3:end) = T2;
T(:,3:3:end) = T3;

s2 = size(T,2);
inds1 = (1:3:s2);
inds2 = (2:3:s2);
inds3 = (3:3:s2);
iskel = [inds1(sk);inds2(sk);inds3(sk)];
iskel = iskel(:);

Tred = T(:,iskel);
tmat = Tred\T;




% so, tmat takes !3! densities on the boundary to a collection 
% of weights and skeleton points which handle farfield scattering
% first component  - 1st density
% second component - 1st density (Hilbert)
% third component  - 2nd density

%%

chnkr.kappa = signed_curvature(chnkr);

[data] = flex_pt(zk,coefs,pxypts,chnkr);
sz = size(data);
T = data.';

T1 = T(:,1:2:end);
T2 = T(:,2:2:end);

tol = 1E-10;
I = eye(size(T1,1));
[sk1,~,~] = id(T1,tol);
[Q1,~] = qr(T1(:,sk1));
T2red = (I-Q1*Q1')*T2;
[sk2,~,~] = id(T2red,tol/norm(T2red,'fro'));
sk_in = unique([sk1,sk2]);

s2 = size(T,2);
inds1 = (1:2:s2);
inds2 = (2:2:s2);
iskel_in = [inds1(sk_in);inds2(sk_in)];
iskel_in = iskel_in(:);
Tred = T(:,iskel_in);
tmat_in = Tred\T;

% so, tmat takes !2! evaluations of a boundary conditions at a
% collection of skeleton points to the boundary conditions on the 
% entire boundary
% first component  - 1st bc
% second component - 2nd bc


%% calculate the reduced scattering matrix ttS

S = inv(Sinv);
nS = size(S,1);

H = chunkermat(chnkr, hilbert, opts2);  


zz = zeros(3*nS/2,nS);
zz(1:3:end,1:2:end) = eye(nS/2);
zz(3:3:end,2:2:end) = eye(nS/2);
zz(2:3:end,1:2:end) = H;
tS = zz*S;
ttS = tmat*tS*tmat_in.';

ttmp = S*tmat_in.';


%% 

% get the reduced bc (on skeleton)
targs = [];
targs.r = chnkr.r(:,sk_in);
targs.d = chnkr.d(:,sk_in);
targs.n = chnkr.n(:,sk_in);
targs.kappa = chnkr.kappa(sk_in);
[data] = flex_pt(zk,coefs,source_loc,targs);

tdens = ttmp*data;
tout = zz*ttmp*data;

v1 = tout(1:3:end);
v2 = tout(2:3:end);
v3 = tout(3:3:end);

targ = [];
targ.r = [0;3];

ch1src = [];
ch1src.r = chnkr.r(:,:);
ch1src.n = chnkr.n(:,:);
ch1src.d = chnkr.d(:,:);
ch1src.w = chnkr.wts(:);
ik1 = ikern1(ch1src,targ)*(v1.*ch1src.w);
ik2 = ikern2(ch1src,targ)*(v3.*ch1src.w);
ik3 = ikern3(ch1src,targ)*(v2.*ch1src.w);

targs = targ.r;
uscat = 0;
count = 0;

rho1 = sol(count+1:2:count+2*chnkr.npt);
rho2 = sol(count+2:2:count+2*chnkr.npt);
count = count + 2*chnkrs(i).npt;

H = chunkermat(chnkr, hilbert, opts2);                                              % Assemble hilbert transforms

uscat =uscat + chunkerkerneval(chnkr, ikern1,rho1, targs) +...
    chunkerkerneval(chnkr, ikern3,H*rho1, targs) + ...
    chunkerkerneval(chnkr, ikern2, rho2, targs);

% calculate reduced outgoing density
tdens = ttS*data;

v1 = tdens(1:3:end);
v2 = tdens(2:3:end);
v3 = tdens(3:3:end);

targ = [];
targ.r = [0;3];

ch1src = [];
ch1src.r = chnkr.r(:,sk);
ch1src.n = chnkr.n(:,sk);
ch1src.d = chnkr.d(:,sk);
ch1src.w = chnkr.wts(sk.');
ik1b = ikern1(ch1src,targ)*(v1.*ch1src.w);
ik2b = ikern2(ch1src,targ)*(v3.*ch1src.w);
ik3b = ikern3(ch1src,targ)*(v2.*ch1src.w);

% check the error
ik1b+ik2b+ik3b-uscat

%% skeleton

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

[sout] = skel_chnkr(chnkr,nu,zk,coefs,kerns,pxypts);

%% now test again

targs = [];
targs.r = sout.sk_in_r;
targs.d = sout.sk_in_d;
targs.n = sout.sk_in_n;
targs.kappa = sout.kappa;
[data] = flex_pt(zk,coefs,source_loc,targs);

tdens2 = sout.scat_red*data;

v1 = tdens2(1:3:end);
v2 = tdens2(2:3:end);
v3 = tdens2(3:3:end);

targ = [];
targ.r = [0;3];

ch1src = [];
ch1src.r = sout.sk_ou_r;
ch1src.n = sout.sk_ou_n;
ch1src.d = sout.sk_ou_d;
ch1src.w = sout.sk_ou_w;
ik1b = ikern1(ch1src,targ)*(v1.*ch1src.w);
ik2b = ikern2(ch1src,targ)*(v3.*ch1src.w);
ik3b = ikern3(ch1src,targ)*(v2.*ch1src.w);

ik1b + ik2b + ik3b-uscat
