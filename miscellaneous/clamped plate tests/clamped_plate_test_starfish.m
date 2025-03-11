clearvars; close all;
clear
iseed = 8675309;
rng(iseed);

% planewave vec

% kvec = 200*[1;-1.5];
% 
% %
% 
% zk = sqrt(norm(kvec));

% discretize domain
zk = 10;
cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 1;
cparams.maxchunklen = 4./zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% plot geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal


tic = start ;
opts = [];

fkern = @(s,t) flex2d.kern(zk, s, t, 'clamped-plate'); 

start = tic;
D = chunkermat(chnkr,fkern, opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)


x = chnkr.d(1,:);
y = chnkr.d(2,:);

denom = sqrt(x.^2+y.^2);
numer = chnkr.d(1,:).*chnkr.d2(2,:)-chnkr.d2(1,:).*chnkr.d(2,:);

kappa = numer./(denom.^3);


A = zeros(2, 2, chnkr.npt);
start = tic;
for i = 1:chnkr.npt
    A(:, :, i) = [-1/2, 0 ; kappa(i), -1/2];
end
t3 = toc; 
fprintf('%5.2e s : time to construct jump matrix\n',t3);                                   % Needs to find a way to aviod this              

K = num2cell(A, [1 2]);
M = blkdiag(K{:}); 
 

[y1, grad, ~, ~, ~] = flex2d.hkdiffgreen(zk, [0;0], chnkr.r);

nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);

normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         

firstbc = 1/(2*zk^2).*y1;
secondbc = 1/(2*zk^2).*normalderiv;  

[nt, ~] = size(D);
lhs = M + D;
rhs = zeros(nt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;



tic;
%sol = gmres(lhs, rhs, [], 1e-13, 200);
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = (-2:0.05:2) + randn()*0.01;                                     % generate some targets
ys = (-2:0.05:2) + randn()*0.01;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);


tic;
in = chunkerinterior(chnkr, targets); 
out = ~in; 
toc



ikern1 = @(s,t) flex2d.kern(zk, s, t, 'first kernel');                              % build the kernel of evaluation          
ikern2 = @(s,t) flex2d.kern(zk, s, t, 'second kernel');

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


true_sol = nan(na, 1);
utarg = nan(na, 1);

[val, ~] = flex2d.hkdiffgreen(zk, [0;0], targets(:,out));

trueval = 1/(2*zk^2).*val;

utarg(out) = Dsol;
true_sol(out) = trueval;

uerr = abs(utarg - true_sol)/max(abs(true_sol(:)));
uerr = reshape(uerr,size(X));

figure(2)
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
title("Absolute error(check-clamped plate kernel on starfish k = 18)", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
%colormap(hsv)
colorbar











