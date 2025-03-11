clear  
clearvars
close all


% kvec = 200*[1;-1.5];
% 
% zk = sqrt(norm(kvec)); 

zk = 0.1;                                                % wave number

cparams= [];
cparams.eps = 1e-6;
%cparams.nover = 0 ;
cparams.maxchunklen = min(0.2, 4./zk);  

           
start = tic;
chnkr = chunkerfunc(@(t) ellipse(t), cparams);            % build the geometry
t1 = toc(start);
fprintf('%5.2e s : time to build geo\n',t1);
                    

figure(1)                                                          % plot the geometry                   
clf
plot(chnkr, '-x')
hold on
quiver(chnkr)
axis equal


opts = [];
opts.sing = 'log';
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped-plate');           % build the desired kernel

start = tic;
D = chunkermat(chnkr,fkern, opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)



kappa = signed_curvature(chnkr);
kappa = kappa(:);

A = zeros(2, 2, chnkr.npt);
start = tic;
for i = 1:chnkr.npt
    A(:, :, i) = [1/2, 0 ; -kappa(i), 1/2];
end
t3 = toc; 
fprintf('%5.2e s : time to construct jump matrix\n',t3);

K = num2cell(A, [1 2]);
M = blkdiag(K{:}); 
 

[y1, grad, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, [3;2], chnkr.r);

[y1K, gradK, ~ , ~ ,~] = chnk.flex2d.helmdiffgreen(zk*(1i), [3;2], chnkr.r);



nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);

normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         
normalderivK = gradK(:, :, 1).*(nx.') + gradK(:, :, 2).*(ny.');

firstbc = 1/(2*zk^2).*y1 - 1/(2*zk^2).*y1K;
secondbc = 1/(2*zk^2).*normalderiv - 1/(2*zk^2).*normalderivK;

[nt, ~] = size(D);
lhs = M + D;

rhs = zeros(nt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;


tic
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = -4:0.05:4;                                    % generate some targets
ys = -4:0.05:4;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = chunkerinterior(chnkr, targets); 
out = ~in; 
toc

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'first kernel');                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'second kernel');


start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets(:, in)) + chunkerkerneval(chnkr, ikern2, rho2, targets(:,in));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

[val, ~] = chnk.flex2d.helmdiffgreen(zk, [3;2], targets(:,in));

[valK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), [3;2], targets(:,in));


trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

utarg(in) = Dsol;
true_sol(in) = trueval;

uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));

figure(2)
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
title("Absolute error ", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colormap(hsv)
colorbar

