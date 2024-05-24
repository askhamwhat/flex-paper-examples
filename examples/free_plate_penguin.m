%clear
%clc

zk = 160;
nu = 1/3;

cparams.eps = 1e-6;
cparams.maxchunklen = 4 ./ zk;       % setting a chunk length helps when the
                                    % frequency is known'
theta = pi/3;
d = -[cos(theta) sin(theta)];


cparams.maxchunklen = 4.0 / zk; % setting a chunk length helps when the
                              % frequency is known
cparams.rounded = 'true';

verts = readmatrix('best_penguin.txt') / 550;

chnkr = chunkerpoly(verts',cparams);
chnkr = refine(chnkr);

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
clf
plot(chnkr, '-x')
title('Chunkr object')
hold on
quiver(chnkr)
axis equal
drawnow


coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel


fkern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
fkern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);   
hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
double = @(s,t) chnk.lap2d.kern(s,t, 'd');

fkern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part', coefs);                     % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

fkern4 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 second part', coefs);                    % kernels in K21 needs to multiply by curvature

fkern5 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);                   % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature

fkern6 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K22 second part', coefs);                    % kernels in K22 needs to multiply by curvature


start = tic;
sysmat = chunkermat(chnkr,fkern, opts);
sysmat1 = chunkermat(chnkr, fkern1, opts);
sysmat2 = chunkermat(chnkr, fkern2, opts);
K21 = chunkermat(chnkr, fkern3, opts);
K21second = chunkermat(chnkr, fkern4, opts);
K21hilbert = chunkermat(chnkr, fkern5, opts);

K22second = chunkermat(chnkr, fkern6, opts);

D = chunkermat(chnkr, double, opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

opts2 = [];
opts2.sing = 'pv';

start = tic;
H = chunkermat(chnkr, hilbert, opts2);                                              % Assemble hilbert transforms
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

kappa = signed_curvature(chnkr);
kappa = kappa(:);
kappa = kappa';

hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);

sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa.').*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + (kappa.').*(K22second);


A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)

M = kron(eye(chnkr.npt), A);

lhs =  M + sysmat ;

return

[~, ~, hess, third] = planewave1(zk, chnkr.r(:,:), d);


nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr);
kappa = kappa(:)';

firstbc = (hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
           coefs(1).*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy));

secondbc = (third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) +...
       third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny))  + ...
        (2-coefs(1)).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:,3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:,4).*(tauy.*tauy.*ny)) + ...
        (1-coefs(1)).*kappa'.*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy+...
        -(hess(:,1).*nx.*nx + hess(:,2).*(2*nx.*ny) + hess(:,3).*ny.*ny));

[nt, ~] = size(sysmat);

rhs = zeros(nt, 1); 
rhs(1:2:end) = -firstbc ; 
rhs(2:2:end) = -secondbc ;

tic
sol = lhs\rhs;
toc;


rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = -6:0.01:6;                                    % generate some targets
ys = -6:0.01:6;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

in = chunkerinterior(chnkr, targets); 
out = ~in; 



ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));


start1 = tic;
uscat = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

[uinc,~, ~] = planewave1(zk,targets(:,out),d);


us = zeros(na, 1);
ui = zeros(na,1);
us(out) = uscat;
ui(out) = uinc;
us = reshape(us,[numel(xs) numel(ys)]);
ui = reshape(ui,[numel(xs) numel(ys)]);
utot = us + ui;


figure(3)
tiledlayout(1,3)
nexttile
h = pcolor(X,Y,real(ui));
h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
title("Incident field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
clim([-1,1]);
colorbar

nexttile
h = pcolor(X,Y,real(us));
h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
title("Scattered field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
clim([-1,1]);
colorbar

nexttile
h = pcolor(X,Y,real(utot));
h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
title("Total field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
clim([-1,1]);
colorbar


