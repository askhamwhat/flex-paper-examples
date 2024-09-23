function [sout] = skel_chnkr(chnkr,nu,zk,coefs,kerns,pxysrf)

sout = [];

fkern = kerns{1};
fkern1= kerns{2};
fkern2= kerns{3};
hilbert = kerns{4};
double  = kerns{5};
fkern3  = kerns{6};
fkern4  = kerns{7};
fkern5  = kerns{8};
fkern6  = kerns{9};
fkern1_trans = kerns{10};
fkern3_trans = kerns{11};

ikern1 = kerns{12};
ikern2 = kerns{13};
ikern3 = kerns{14};

% build the scattering matrix

src_chnkr = chnkr;

opts = [];
opts.sing = 'log';

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

S = inv(Sinv);
sout.S = S;



%%%%%%%%% outgoing skeleton

pt_mins = min(chnkr.r(:,:),[],2);
pt_maxs = max(chnkr.r(:,:),[],2);

cent = (pt_mins+pt_maxs)/2;
rad  = max(pt_maxs-pt_mins)/2;

pxypts = cent+(pxysrf*rad);
%pxypts = pxysrf;
T1 = chunkerkernevalmat(chnkr, ikern1, pxypts);
T2 = chunkerkernevalmat(chnkr, ikern3, pxypts);
T3 = chunkerkernevalmat(chnkr, ikern2, pxypts);

T = [T1;T2;T3];
tol = 1E-10;
I = eye(size(T1,1));
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

sout.tmat  = tmat;
sout.sk    = sk;
sout.iskel = iskel;

sout.sk_ou_r = chnkr.r(:,sk);
sout.sk_ou_d = chnkr.d(:,sk);
sout.sk_ou_n = chnkr.n(:,sk);
sout.sk_ou_w = chnkr.wts(sk.');

%%%%%%%% incoming skeleton

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

sout.Tmat_in = tmat_in;
sout.sk_in = sk_in;
sout.iskel_in = iskel_in;
sout.kappa = chnkr.kappa(sk_in);
sout.sk_in_r = chnkr.r(:,sk_in);
sout.sk_in_d = chnkr.d(:,sk_in);
sout.sk_in_n = chnkr.n(:,sk_in);

%%%%%%% post processing

nS = size(S,1);

H = chunkermat(chnkr, hilbert, opts2);  


zz = zeros(3*nS/2,nS);
zz(1:3:end,1:2:end) = eye(nS/2);
zz(3:3:end,2:2:end) = eye(nS/2);
zz(2:3:end,1:2:end) = H;
tS = zz*S;
ttS = tmat*tS*tmat_in.';

sout.scat_red = ttS;
sout.expand   = S*tmat_in.';
sout.pxypts = pxypts;
end