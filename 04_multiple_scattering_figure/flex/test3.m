% Parameters

n_objs = 2; % number of objects
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

%%
%%%%% build the full translation matrix

src_chnkr = chnkrs{1};
targ_chnkr= chnkrs{2};
wts = src_chnkr.wts(:)';

sysmat = zeros(targ_chnkr.npt*2,src_chnkr.npt*3);

fsys = fkern(src_chnkr,targ_chnkr);
sysmat(:,1:3:end) = fsys(:,1:2:end).*wts;
sysmat(:,3:3:end) = fsys(:,2:2:end).*wts;

sysmat1 = fkern1_trans(src_chnkr,targ_chnkr).*wts;
sysmat2 = fkern2(src_chnkr,targ_chnkr).*wts;

K21 = fkern3_trans(src_chnkr,targ_chnkr).*wts;
K21second = fkern4(src_chnkr,targ_chnkr).*wts;
K21hilbert = fkern5(src_chnkr,targ_chnkr).*wts;
K22second = fkern6(src_chnkr,targ_chnkr).*wts;


kappa = signed_curvature(targ_chnkr);
kappa = kappa(:)';

hilb = sysmat1;
hilb2 = sysmat2;

mat1 = sysmat(1:2:end, 1:3:end);
mat4 = sysmat(2:2:end, 3:3:end);

sysmat(1:2:end, 1:3:end) = mat1;
sysmat(2:2:end, 1:3:end) = K21  + (kappa.').*(K21second);
sysmat(2:2:end, 3:3:end) = mat4 + (kappa.').*(K22second);

sysmat(1:2:end,2:3:end) = hilb;
sysmat(2:2:end,2:3:end) = hilb2 + (kappa.').*(K21hilbert);

T =  sysmat;

%%
%%% checking the skeleton, expand, and eval matrices

iskel_ou = sout{1}.iskel;
iskel_in = sout{2}.iskel_in;
texpand = sout{1}.tmat;
tevalua = sout{2}.Tmat_in;
Tred = T(iskel_in,iskel_ou);
Tfull = tevalua.'*Tred*texpand;

%%
%%% now checking the fast Tred generation...
src = [];
src.r = sout{1}.sk_ou_r;
src.n = sout{1}.sk_ou_n;
src.d = sout{1}.sk_ou_d;
src.wts = sout{1}.sk_ou_w;

tar = []; 
tar.r = sout{2}.sk_in_r;
tar.n = sout{2}.sk_in_n;
tar.d = sout{2}.sk_in_d;
tar.kappa = sout{2}.kappa;

[sysmat] = get_trans_mat(src,tar,kerns);
Tfull = tevalua.'*sysmat*texpand;
norm(T-Tfull,'fro')
return