zk = 1;               % our k (wave number)
nu = 1/3;
cparams = [];
cparams.maxchunklen = 2;

chnkr1 = chunkerfunc(@(t) ellipse(t), cparams);
chnkr2 = chunkerfunc(@(t) ellipse(t), cparams);
center = [4; 4];
position = [4; 4];
chnkr2 = chnkr2.move([0;0],center,0,1);

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
tiledlayout(1,3,'TileSpacing','tight')
nexttile
plot(chnkr1, '-x')
hold on
plot(chnkr2, '-x')
title('Chunkr objects')
hold on
quiver(chnkr1); hold on
quiver(chnkr2)
axis equal
xlim([-4 8])
ylim([-4 8])
drawnow

coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate first part', coefs);                          % build the desired kernel
fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
fkern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);                    % first part K21 coupled with hilbert transforms

hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
double = @(s,t) chnk.lap2d.kern(s, t, 'd');

fkern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part', coefs);             % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)
fkern4 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 second part', coefs);            % kernels in K21 needs to multiply by curvature
fkern5 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);           % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature
fkern6 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K22 second part', coefs);            % kernels in K22 needs to multiply by curvature

start = tic;

% Scattering matrix S1 

sysmat = chunkermat(chnkr1,fkern, opts);
sysmat1 = chunkermat(chnkr1, fkern1, opts);
sysmat2 = chunkermat(chnkr1, fkern2, opts);

mat3 = chunkermat(chnkr1, fkern3, opts);
K21second = chunkermat(chnkr1, fkern4, opts);
K21hilbert = chunkermat(chnkr1, fkern5, opts);
K22second = chunkermat(chnkr1, fkern6, opts);

D = chunkermat(chnkr1, double, opts);

opts2 = [];
opts2.sing = 'pv';

H = chunkermat(chnkr1, hilbert, opts2);                                              % Assemble hilbert transforms

kappa = signed_curvature(chnkr1);
kappa = kappa(:)';

hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = mat3 +  hilb2 + diag(kappa)*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + diag(kappa)*(K22second);

A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
M = kron(eye(chnkr1.npt), A);
S1inv =  M + sysmat ;

% Scattering matrix S2

sysmat = chunkermat(chnkr2,fkern, opts);
sysmat1 = chunkermat(chnkr2, fkern1, opts);
sysmat2 = chunkermat(chnkr2, fkern2, opts);

mat3 = chunkermat(chnkr2, fkern3, opts);
K21second = chunkermat(chnkr2, fkern4, opts);
K21hilbert = chunkermat(chnkr2, fkern5, opts);
K22second = chunkermat(chnkr2, fkern6, opts);

D = chunkermat(chnkr2, double, opts);
H = chunkermat(chnkr2, hilbert, opts2);                                              % Assemble hilbert transforms

kappa = signed_curvature(chnkr2);
kappa = kappa(:)';

hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = mat3 +  hilb2 + diag(kappa)*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + diag(kappa)*(K22second);

A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
M = kron(eye(chnkr2.npt), A);
S2inv =  M + sysmat;


% Translation matrix T12

fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert unsubtract', coefs);  % K11 without the hilbert subtraction
fkern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part unsubtract', coefs);             % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

wts = chnkr2.wts(:)';

sysmat = fkern(chnkr2,chnkr1);
sysmat(:,1:2:end) = sysmat(:,1:2:end)*diag(wts);
sysmat(:,2:2:end) = sysmat(:,2:2:end)*diag(wts);
sysmat1 = fkern1(chnkr2,chnkr1)*diag(wts);
sysmat2 = fkern2(chnkr2,chnkr1)*diag(wts);

mat3 = fkern3(chnkr2,chnkr1)*diag(wts);
K21second = fkern4(chnkr2,chnkr1)*diag(wts);
K21hilbert = fkern5(chnkr2,chnkr1)*diag(wts);
K22second = fkern6(chnkr2,chnkr1)*diag(wts);

H = chunkermat(chnkr2, hilbert, opts2);                                              % Assemble hilbert transforms

kappa = signed_curvature(chnkr1);
kappa = kappa(:)';

hilb = sysmat1*H ;
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = mat3 +  hilb2 + diag(kappa)*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + diag(kappa)*(K22second);

T12 =  sysmat;

% Translation matrix T21

wts = chnkr1.wts(:)';

sysmat = fkern(chnkr1,chnkr2);
sysmat(:,1:2:end) = sysmat(:,1:2:end)*diag(wts);
sysmat(:,2:2:end) = sysmat(:,2:2:end)*diag(wts);
sysmat1 = fkern1(chnkr1,chnkr2)*diag(wts);
sysmat2 = fkern2(chnkr1,chnkr2)*diag(wts);

mat3 = fkern3(chnkr1,chnkr2)*diag(wts);
K21second = fkern4(chnkr1,chnkr2)*diag(wts);
K21hilbert = fkern5(chnkr1,chnkr2)*diag(wts);
K22second = fkern6(chnkr1,chnkr2)*diag(wts);

D = chunkermat(chnkr1, double, opts);
H = chunkermat(chnkr1, hilbert, opts2);                                              % Assemble hilbert transforms

kappa = signed_curvature(chnkr2);
kappa = kappa(:)';

hilb = sysmat1*H ;
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = mat3 +  hilb2 + diag(kappa)*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + diag(kappa)*(K22second);

T21 =  sysmat;

% Assembling block matrix

lhs = [S1inv T12; T21 S2inv];

t1 = toc(start);
fprintf('%5.2f s : time to assemble matrix\n',t1)

% Evaluating righthand side(s)

start = tic;

[~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, position, chnkr1.r);
[~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), position, chnkr1.r);

nx = chnkr1.n(1,:).'; 
ny = chnkr1.n(2,:).';

dx = chnkr1.d(1,:).';
dy = chnkr1.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr1);
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

rhs = zeros(2*chnkr1.npt+2*chnkr2.npt, 1); 
rhs(1:2:2*chnkr1.npt) = firstbc ; 
rhs(2:2:2*chnkr1.npt) = secondbc;

[~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, position, chnkr2.r);
[~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), position, chnkr2.r);

nx = chnkr2.n(1,:).'; 
ny = chnkr2.n(2,:).';

dx = chnkr2.d(1,:).';
dy = chnkr2.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr2);
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

rhs(2*chnkr1.npt+1:2:end) = firstbc ; 
rhs(2*chnkr1.npt+2:2:end) = secondbc ;

t1 = toc(start);
fprintf('%5.2f s : time to evaluate BCs \n',t1) 

% correct_density = S1inv \ rhs(1:2*chnkr1.npt);

sol = lhs\rhs;

rho1 = sol(1:2:2*chnkr1.npt); 
rho2 = sol(2:2:2*chnkr1.npt); 
rho3 = sol(2*chnkr1.npt+1:2:end); 
rho4 = sol(2*chnkr1.npt+2:2:end);

% correct_rho1 = correct_density(1:2:end); 
% correct_rho2 = correct_density(2:2:end); 

xs = -4:0.2:8;                                    % generate some targets
ys = -4:0.2:8;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = chunkerinterior(chnkr1, targets)+chunkerinterior(chnkr2, targets); 
out = ~in; 
toc


ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

H = chunkermat(chnkr1, hilbert, opts2);                                              % Assemble hilbert transforms

start1 = tic;
Dsol1 = chunkerkerneval(chnkr1, ikern1,rho1, targets(:, out)) +...
    chunkerkerneval(chnkr1, ikern3,H*rho1, targets(:, out)) + ...
    chunkerkerneval(chnkr1, ikern2, rho2, targets(:,out));

H = chunkermat(chnkr2, hilbert, opts2);                                              % Assemble hilbert transforms

Dsol2 = chunkerkerneval(chnkr2, ikern1,rho3, targets(:, out)) + ...
    chunkerkerneval(chnkr2, ikern3,H*rho3, targets(:, out)) + ...
    chunkerkerneval(chnkr2, ikern2, rho4, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


[val, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, position, targets(:,out));
[valK, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), position, targets(:,out));

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

true_sol = zeros(na, 1);
utarg = zeros(na, 1);
% correct_utarg = zeros(na,1);

utarg(out) = Dsol1 + Dsol2;
true_sol(out) = trueval;

% Dsol1 = chunkerkerneval(chnkr1, ikern1,correct_rho1, targets(:, out)) +...
%     chunkerkerneval(chnkr1, ikern3,H*correct_rho1, targets(:, out)) + ...
%     chunkerkerneval(chnkr1, ikern2, correct_rho2, targets(:,out));
% 
% correct_utarg(out) = Dsol1;

utarg = reshape(utarg,size(X));
% correct_utarg = reshape(correct_utarg,size(X));
true_sol = reshape(true_sol,size(X));

% nexttile
% h = pcolor(X,Y,real(utarg));
% set(h,'EdgeColor','None'); hold on;
% title("Coupled BIE solution")
% plot(chnkr1,'w-','LineWidth',2);
% hold on
% plot(chnkr2,'w-','LineWidth',2);
% colorbar
% 
% h = pcolor(X,Y,real(correct_utarg));
% set(h,'EdgeColor','None'); hold on;
% title("Independent BIE solution")
% plot(chnkr1,'w-','LineWidth',2);
% hold on
% plot(chnkr2,'w-','LineWidth',2);
% colorbar

nexttile
h = pcolor(X,Y,real(utarg));
set(h,'EdgeColor','None'); hold on;
title("True solution")
plot(chnkr1,'w-','LineWidth',2);
hold on
plot(chnkr2,'w-','LineWidth',2);
colorbar

nexttile
uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));
h = pcolor(X,Y,log10(abs(uerr) / max(abs(true_sol),[],'all')));
set(h,'EdgeColor','None'); hold on;
title("Log_{10} relative error (coupled)")
plot(chnkr1,'w-','LineWidth',2);
hold on
plot(chnkr2,'w-','LineWidth',2);
colorbar

return

% nexttile
% uerr = correct_utarg - true_sol;
% uerr = reshape(uerr,size(X));
% h = pcolor(X,Y,log10(abs(uerr) / max(abs(true_sol),[],'all')));
% set(h,'EdgeColor','None'); hold on;
% title("Log_{10} relative error (individual)")
% plot(chnkr1,'w-','LineWidth',2);
% hold on
% plot(chnkr2,'w-','LineWidth',2);
% colorbar
 

% Debugging section

K11 = mat1 + hilb;
K12 = sysmat(1:2:end,2:2:end);
K21 = mat3 +  hilb2 + diag(kappa)*(K21hilbert*H + K21second);
K22 = mat4 + diag(kappa)*(K22second);

testbc1 = K11*correct_rho1 + K12*correct_rho2;
testbc2 = K21*correct_rho1 + K22*correct_rho2;

figure(2)
tiledlayout(1,2)
nexttile
title('First BC residual error')
plot(real(testbc1 - firstbc),'o')

nexttile
title('Second BC residual error')
plot(real(testbc2 - secondbc ),'o')


% title('First density error')
% plot(real(correct_rho1),'o')
% 
% nexttile
% title('Second BC residual error')
% plot(real(correct_rho2 ),'o')

