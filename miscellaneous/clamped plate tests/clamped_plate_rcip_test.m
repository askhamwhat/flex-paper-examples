clear 
clc
verts = [0, 3,  0, 3; 
         0, 0,  3, 3;];

% kvec = 100*[1;-1.5];

zk = 18;   

cparams= [];
cparams.eps = 1e-10;
cparams.maxchunklen = 4.0./zk;

edge2verts = [ -1 1 0 0;                    %  edge2verts (connect the geometry with the line)
               1  0 -1 0;
               0  0  1 -1;
               0  -1  0  1]; 

fchnks{1} = @(t) rotateline(t, 0);
fchnks{2} = @(t) rotateline(t, pi/2);
fchnks{3} = @(t) rotateline(t, pi);
fchnks{4} = @(t) rotateline(t, 3*pi/2);


cgrph = chunkgraph(verts, edge2verts, fchnks, cparams);

% 
% plot(cgrph); hold on; 
% quiver(cgrph);

opts = [];
opts.rcip = true;

opts.nsub_or_tol = 40;
fkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped-plate-rcip');           % build the desired kernel
start = tic;
D = chunkermat(cgrph,fkern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)




[y1, grad, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, [2;2], chnkr.r);

[y1K, gradK, ~ , ~ ,~] = chnk.flex2d.helmdiffgreen(zk*(1i), [2;2], chnkr.r);



nx = cgrph.n(1,:); 
ny = cgrph.n(2,:);

normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC                          
normalderivK = gradK(:, :, 1).*(nx.') + gradK(:, :, 2).*(ny.');

firstbc = 1/(2*zk^2).*y1 - 1/(2*zk^2).*y1K;
secondbc = 1/(2*zk^2).*normalderiv - 1/(2*zk^2).*normalderivK;                              % Dirichlet and Neumann BC                          



M = eye(2*cgrph.npt);

[nt, ~] = size(D);

lhs = M + D;
rhs = zeros(nt, 1); 
rhs(1:2:end) = -2.*y1 ; 
rhs(2:2:end) = -2.*normalderiv;


tic
%sol = gmres(lhs, rhs, [], 1e-13, 400);
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        

zz = length(rho1);

rho11 = rho1(1:zz/4);                                   % first density ( cgrph.echnks(1))
rho12 = rho1(zz/4 + 1:zz/2);                            % first density ((cgrph.echnks(2))
rho13 = rho1(zz/2 +1 : 3*zz/4);                         % first density ((cgrph.echnks(3))
rho14 = rho1(3*zz/4 +1 :end);                           % first density ((cgrph.echnks(4))

rho21 = rho2(1:zz/4);                                    % second density ( cgrph.echnks(1))
rho22 = rho2(zz/4 +1 : zz/2);                           % second density ( cgrph.echnks(2))
rho23 = rho2(zz/2 +1 : 3*zz/4);                         % second density ( cgrph.echnks(3))
rho24 = rho2(3*zz/4 +1 :end);                           % second density ( cgrph.echnks(4))




xs = -4:0.1:4;                                    % generate some targets
ys = -4:0.1:4;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
in = ((X >= 0 & X<= 3) & (Y >= 0 & Y<= 3)) ; 
out = ~in; 
toc

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'first kernel');                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'second kernel');

start1 = tic;
Dsol = chunkerkerneval(cgrph.echnks(1), ikern1,rho11, targets(:, out)) +...
    chunkerkerneval(cgrph.echnks(1), ikern2, rho21, targets(:, out)) + ...
    chunkerkerneval(cgrph.echnks(2), ikern1, rho12, targets(:, out)) + ...
    chunkerkerneval(cgrph.echnks(2), ikern2, rho22, targets(:, out)) + ...
    chunkerkerneval(cgrph.echnks(3), ikern1, rho13, targets(:, out)) + ...
    chunkerkerneval(cgrph.echnks(3), ikern2, rho23, targets(:, out))+ ...
    chunkerkerneval(cgrph.echnks(4), ikern1, rho14, targets(:, out)) + ...
    chunkerkerneval(cgrph.echnks(4), ikern2, rho24, targets(:, out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t1)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

[val, ~] = chnk.flex2d.helmdiffgreen(zk, [2;2], targets(:,out));

[valK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), [2;2], targets(:,out));


trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;


utarg(out) = Dsol;
true_sol(out) = trueval;

uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));

figure(2)
h = pcolor(X,Y,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
title("Absolute error (check-clamped plate kernel on a square) using rcip", 'FontSize',16)
plot(cgrph,'w-','LineWidth',2);
axis equal
colorbar








