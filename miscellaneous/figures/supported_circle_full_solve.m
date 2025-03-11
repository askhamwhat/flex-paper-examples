clear 
clearvars

zk = 5;
nu = 0.3;

cparams = [];
pref = [];
pref.k = 16;
cparams.eps = 1e-5;
cparams.nover = 0;
cparams.maxchunklen = 0.25; % 4./zk;       % setting a chunk length helps when the
                                    % frequency is known'
R = 1;
chnkr = chunkerfunc(@(t) circle(t), cparams, pref);
chnkr = chnkr.move([0;0],[0;0],0,R);
chnkr = chnkr.sort();
centre = [0.3; 0];
coefs = nu;

kappa = signed_curvature(chnkr);
kappa = kappa(:);

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
clf
plot(chnkr, '-x')
hold on
quiver(chnkr)
hold on 
%scatter(xs(:),ys(:),36,kp,'filled')
axis equal
drawnow

ikern1 =  @(s,t) suppkern(zk, s, t, 'supported plate',coefs);           % build the desired kernel
ikern2 =  @(s,t) suppkern(zk, s, t, 'supported plate K21',coefs);           % build the desired kernel

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

start = tic;
M = chunkermat(chnkr,ikern1, opts);
M2 = chunkermat(chnkr,ikern2, opts2);
M2(isnan(M2)) = 0;

K11fo = M(1:2:end,1:2:end);
K12fo = M(1:2:end,2:2:end);
K21fo = M(2:2:end,1:2:end);
K22fo = M(2:2:end,2:2:end);

% Ellipse:  

c = 1;

xs = chnkr.r(1,:,:);
ys = chnkr.r(2,:,:);

xs = xs(:);
ys = ys(:);

t = atan2(c*ys,xs);

x1 = -c*sin(t);    
y1 = cos(t);

x2 = -c*cos(t);
y2 = -sin(t);

x3 = c*sin(t);
y3 = -cos(t);

x4 = c*cos(t);
y4 = sin(t);


kp2 = (-x4.*y1.^3 + y1.^2.*(4*y3.*x2+y4.*x1) + x1.*(-3*x2.^2.*y2 + y4.*x1.^2 - 4*x3.*x1.*y2 - 3*y2.^3) ...
     + y1.*(3*x2.^3 - x1.*(x4.*x1 + 4*y3.*y2) + x2.*(4*x3.*x1 + 3*y2.^2))) ./ (x1.^2 + y1.^2).^(7/2) ...
     - 6*(x1.*x2 + y1.*y2).*((x1.^2+y1.^2).*(x1.*y3 - y1.*x3) + 3*(y1.*x2 - x1.*y2).*(x1.*x2 + y1.*y2)) ./ (x1.^2 + y1.^2).^(9/2);

k21diag = (nu - 1)*(12*kappa.^3*(nu^2 - nu + 4) + kp2*(-5*nu^2 + 4*nu + 33))/(48*pi*(nu - 3)) + 1i*zk^2/64*(3+nu)*((-1+nu)*(7+nu)/(3-nu))*kappa;

K21fo = M2 + diag(k21diag.*chnkr.wts(:));

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

[val, ~, hess, ~, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);

firstbc = 1/(2*zk^2).*val ;

secondbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
           nu/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy));

[nt, ~] = size(lhs);

rhs = zeros(nt, 1); 
rhs(1:2:end) = firstbc ; 
rhs(2:2:end) = secondbc;

start = tic;
sol = lhs\rhs;
t3 = toc(start); 
fprintf('%5.2e s : time to solve system \n',t3);



xs = centre(1)-4-0.05:0.2:centre(1)+4+0.05;                                    % generate some targets
ys = centre(2)-4-0.05:0.2:centre(2)+4+0.05;  
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = chunkerinterior(chnkr, targets); 
out = ~in; 
toc

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);  % second density


ikern1 = @(s,t) suppkern(zk, s, t, 'supported plate K1 eval',coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) suppkern(zk, s, t, 'supported plate K2 eval',coefs);

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1, rho1, targets(:, out)) + ...
        chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);


[val, ~] = flex2d.hkdiffgreen(zk, centre, targets(:,out));

trueval = 1/(2*zk^2).*val ;



utarg(out) = Dsol;
true_sol(out) = trueval;
utarg = reshape(utarg,size(X));
true_sol = reshape(true_sol,size(X));
uerr = utarg - true_sol;

figure(3);
tiledlayout(1,3);
nexttile
h = pcolor(X,Y,real(true_sol));
set(h,'EdgeColor','None'); hold on;
title("True solution", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile
h = pcolor(X,Y,real(utarg));
set(h,'EdgeColor','None'); hold on;
title("BIE solution", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile
h = pcolor(X,Y,log10(abs(true_sol - utarg)));
set(h,'EdgeColor','None'); hold on;
% set(h,'FaceColor','interp');
title("Absolute Error", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
%colormap(hsv)
colorbar

disp('Maximum relative error:')
disp(max(abs(uerr(:)) / max(abs(true_sol),[],'all')))