clear
clc

% CREATING GEOMETRY  - remember to change kappa prime calculation in kern.m

zk = 3;
nu = 0.3;
cparams = [];
cparams.maxchunklen = 2;       % setting a chunk length helps when the
                                    % frequency is known'
theta = 0;
d = -[cos(theta) sin(theta)];

chnkr = chunkerfunc(@(t) starfish3arms(t), cparams);
chnkr = chnkr.move([0;0],[0;0],pi/3,1);
chnkr = chnkr.sort();
centre = [0; 0];
coefs = nu;

% PLOTTING GEOMETRY

figure(1)                                                   % plot the chunker-object (supposed to be a starfish3arms centered at 1 with radius 1)
clf
tiledlayout(1,3,"Padding","loose", 'TileSpacing','tight')
nexttile(1)
% polarplot(ones(1,200), 'k' )
% set(gca, 'ThetaTick', (0:pi/2:3*pi/2)*180/(pi))
% set(gca, 'ThetaTickLabel', {'0','\pi/2','\pi','3\pi/2'})
plot(chnkr, '-k')
xlim([-3 3])
ylim([-3 3])
title('Scattering object')

axis off

C = orderedcolors("gem");
blue = C(1,:);
red = C(2,:);
orange = C(3,:);

% CLAMPED PLATE

coefs = [nu; 0];
opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped-plate', coefs);        % build the desired kernel


start = tic;
D = chunkermat(chnkr,fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble chnkermat\n',t1)


start = tic;
A = zeros(2, 2, chnkr.npt);
kappa = signed_curvature(chnkr);
kappa = kappa(:);
for i = 1:chnkr.npt
    A(:, :, i) = [-0.5, 0 ; kappa(i), -0.5];
end
t3 = toc; 
fprintf('%5.2e s : time to construct jump matrix\n',t3);

K = num2cell(A, [1 2]);
M = blkdiag(K{:}); 
 
[nt, ~] = size(D);
lhs = M + D;

%[val, grad, hess, third, ~] = chnk.flex2d.hkdiffgreen(1i*zk, [0; 0], chnkr.r);
[val, grad, hess, third] = planewave1(zk, chnkr.r(:,:), d);


nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr);
kappa = kappa(:)';

firstbc = val;
secondbc = grad(:,1).*nx + grad(:,2).*ny;

rhs = zeros(nt, 1); 
rhs(1:2:end) = -firstbc ; 
rhs(2:2:end) = -secondbc ;

tic
sol = lhs\rhs;
toc;


rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        % second density


R = 1000;
thetas = -pi:0.05:pi;
xs = R*cos(thetas);
ys = R*sin(thetas);
targets = [xs; ys];
[~,na] = size(targets);

in = chunkerinterior(chnkr, targets); 
out = ~in; 

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'first kernel');                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'second kernel');

start1 = tic;
ufar = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
ufar = sqrt(R)*ufar*exp(-1i*zk*R);

nexttile(2)
plot(thetas, abs(ufar), 'Color', blue)
xlim([-pi pi])
set(gca, 'XTick', -pi:pi/2:pi)
set(gca, 'XTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'})
hold on 

ylabel('|u^{(s)}|')
title('Magnitude of scattered field')

nexttile(3)
phase = atan2(imag(ufar),real(ufar));
phase(logical((phase > 0).*(thetas < -pi/2)')) = phase(logical((phase > 0).*(thetas < -pi/2)')) - 2*pi;
phase(logical((phase > 0).*(thetas > pi/2)')) = phase(logical((phase > 0).*(thetas > pi/2)')) - 2*pi;
plot(thetas, phase, 'Color', blue)
%plot(thetas, real(ufar)./abs(ufar))
xlim([-pi pi])
set(gca, 'XTick', -pi:pi/2:pi)
set(gca, 'XTickLabel', {'-\pi','-\pi/2','0','\pi/2','\pi'})
set(gca, 'YTick', -3*pi/2:pi/2:pi)
set(gca, 'YTickLabel', {'-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi'})
hold on

ylabel('\phi')
title('Phase of scattered field')


% FREE PLATE


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

thetas = -pi:0.05:pi;
xs = R*cos(thetas);
ys = R*sin(thetas);
targets = [xs; ys];
[~,na] = size(targets);

in = chunkerinterior(chnkr, targets); 
out = ~in; 

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));

start1 = tic;
ufar = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)
ufar = sqrt(R)*ufar*exp(-1i*zk*R);

nexttile(2)
plot(thetas, abs(ufar), 'Color', red)

nexttile(3)
phase = atan2(imag(ufar),real(ufar));
phase((phase > 0)) = phase((phase > 0)) - 2*pi;
plot(thetas, phase, 'Color', red)
%plot(thetas, real(ufar)./abs(ufar))

% SUPPORTED PLATE

coefs = [nu; 0];
ikern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported plate 3',coefs);           % build the desired kernel

opts = [];
opts.sing = 'log';

start = tic;
M = chunkermat(chnkr,ikern1, opts);

K11fo = M(1:2:end,1:2:end);
K12fo = M(1:2:end,2:2:end);
K21fo = M(2:2:end,1:2:end);
K22fo = M(2:2:end,2:2:end);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(1:2:end,1:2:end) = K11fo - 0.5*eye(chnkr.npt) ;
M(2:2:end,1:2:end) = K21fo + c0.*kappa.^2.*eye(chnkr.npt);
M(1:2:end,2:2:end) = K12fo;
M(2:2:end,2:2:end) = K22fo - 0.5*eye(chnkr.npt);

lhs = M;

t3 = toc(start); 
fprintf('%5.2e s : time to assemble lhs \n',t3);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';
ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);



firstbc = val;
secondbc = (hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
           coefs(1).*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy));

[nt, ~] = size(lhs);

rhs = zeros(nt, 1); 
rhs(1:2:end) = -firstbc ; 
rhs(2:2:end) = -secondbc;

start = tic;
sol = lhs\rhs;
t3 = toc(start); 
fprintf('%5.2e s : time to solve system \n',t3);


thetas = -pi:0.05:pi;
xs = R*cos(thetas);
ys = R*sin(thetas);
targets = [xs; ys];
[~,na] = size(targets);


rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);  % second density


ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported plate K1 eval 3',coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported plate K2 eval 3',coefs);

start1 = tic;
ufar = chunkerkerneval(chnkr, ikern1, rho1, targets) + ...
        chunkerkerneval(chnkr, ikern2, rho2, targets);
t2 = toc(start1);

ufar = sqrt(R)*ufar*exp(-1i*zk*R);

nexttile(2)
plot(thetas, abs(ufar), 'Color', orange)
hold on 

nexttile(3)
phase = atan2(imag(ufar),real(ufar));
phase(logical((phase > 0).*(thetas < -pi/2)')) = phase(logical((phase > 0).*(thetas < -pi/2)')) - 2*pi;
phase(logical((phase > 0).*(thetas > pi/2)')) = phase(logical((phase > 0).*(thetas > pi/2)')) - 2*pi;
plot(thetas, phase, 'Color', orange)
%plot(thetas, real(ufar)./abs(ufar))
hold on 



% HELMHOLTZ DIRICHLET


fkern = kernel('helm','c',zk,[1,-zk*1i]);
start = tic; sysmat = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

[nt, ~] = size(sys);

rhs = zeros(nt, 1); 
firstbc = val;
rhs(1:1:end) = -firstbc ; 


sol = sys \ rhs;
fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot
% xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot); 
% ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
% [xxtarg,yytarg] = meshgrid(xtarg,ytarg);
% targets = zeros(2,length(xxtarg(:)));
% targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:)



fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

ufar = chunkerkerneval(chnkr,fkern,sol,targets); t1 = toc(start);
ufar = sqrt(R)*ufar*exp(-1i*zk*R);

nexttile(2)
plot(thetas, abs(ufar),  '--', 'Color', blue)

nexttile(3)
phase = atan2(imag(ufar),real(ufar));
phase(logical((phase > 0).*(thetas < -pi/2)')) = phase(logical((phase > 0).*(thetas < -pi/2)')) - 2*pi;
phase(logical((phase > 0).*(thetas > pi/2)')) = phase(logical((phase > 0).*(thetas > pi/2)')) - 2*pi;
plot(thetas, phase, '--', 'Color', blue)
%plot(thetas, real(ufar) ./ abs(ufar), '--')

% HELMHOLTZ NEUMANN

fkern = @(s,t) chnk.helm2d.kern(zk, s, t,'sprime');
sysmat = chunkermat(chnkr,fkern);


lhs = -0.5*eye(chnkr.k*chnkr.nch) + sysmat;


nx = chnkr.n(1,:).' ; ny = chnkr.n(2,:).';

firstbc = grad(:,1).*nx + grad(:,2).*ny;

rhs = -firstbc;

sol = lhs\rhs;

% compute layer potential based on oversample boundary

fkern2 = @(s,t) chnk.helm2d.kern(zk, s, t,'s');
ufar = chunkerkerneval(chnkr,fkern2,sol,targets); 
ufar = sqrt(R)*ufar*exp(-1i*zk*R);

nexttile(2)
plot(thetas, abs(ufar),  '--','Color', red)

nexttile(3)
phase = atan2(imag(ufar),real(ufar));
phase((phase > 0)) = phase((phase > 0)) - 2*pi;
plot(thetas, phase, '--', 'Color', red)
% plot(thetas, real(ufar) ./ abs(ufar), '--')




legend('Clamped','Free','Supported','Dirichlet (Helmholtz)', 'Neumann (Helmholtz)', 'Location','eastoutside')





%
% 
% maxin = max(abs(uin(:)));
% maxsc = max(abs(uin(:)));
% maxtot = max(abs(uin(:)));
% 
% maxu = max(max(maxin,maxsc),maxtot);
% 
% figure(2)
% clf
% 
% t = tiledlayout(1,3,'TileSpacing','compact');
% 
% nexttile
% zztarg = nan(size(xxtarg));
% zztarg(out) = uin;
% h=pcolor(xxtarg,yytarg,imag(zztarg));
% set(h,'EdgeColor','none')
% clim([-maxu,maxu])
% %colormap(brewermap([],'RdBu'));
% hold on
% plot(chnkr,'k','LineWidth',2)
% axis equal tight
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% title('$u^{\textrm{inc}}$','Interpreter','latex','FontSize',12)
% 
% nexttile
% zztarg = nan(size(xxtarg));
% zztarg(out) = uscat;
% h=pcolor(xxtarg,yytarg,imag(zztarg));
% set(h,'EdgeColor','none')
% clim([-maxu,maxu])
% %colormap(brewermap([],'RdBu'));
% hold on
% plot(chnkr,'k','LineWidth',2)
% axis equal tight
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% 
% title('$u^{\textrm{scat}}$','Interpreter','latex','FontSize',12)
% 
% nexttile
% zztarg = nan(size(xxtarg));
% zztarg(out) = utot;
% h=pcolor(xxtarg,yytarg,imag(zztarg));
% set(h,'EdgeColor','none')
% clim([-maxu,maxu])
% %colormap(brewermap([],'RdBu'));
% hold on
% plot(chnkr,'k','LineWidth',2)
% axis equal tight
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% 
% title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)


% 
% %val = besselh(0,zk*sqrt(targets(1,:).^2 + targets(2,:).^2));
% val = chnk.flex2d.hkdiffgreen(zk, [0; 0], targets) / (2*zk^2);
% val = 4*(2*zk^2)*sqrt(zk*R*pi/2)*val*exp(-1i*zk*R)/1i;
% val = real(val);
% 
% 
% plot(thetas, val)
% ylabel('|u_s|')
% xlabel('\theta')
% xlim([0 2*pi])
% set(gca, 'XTick', 0:pi/2:2*pi)
% set(gca, 'XTickLabel', {'0','\pi/2','\pi','3\pi/2','2\pi'})
% 
% title('Point source (normalization)')



% plot absolute value and phase - phase is real part divided by absolute
% value 