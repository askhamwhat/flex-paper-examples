close all
 
clearvars



% kvec = 200*[1;-1.5];
% 
% %
% 
% zk = sqrt(norm(kvec));               % our k (wave number)
zk = 0.01;
cparams = [];

cparams.eps = 1e-12;
cparams.maxchunklen = min(0.5,4./zk);       % setting a chunk length helps when the
                                    % frequency is known'

chnkr = chunkerfunc(@(t) circle(t), cparams);

figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
clf
plot(chnkr, '-x')
hold on
quiver(chnkr)
axis equal


opts = [];
opts.sing = 'log';
fkern =  @(s,t) flex2d.kern(zk, s, t, 'clamped-plate');           % build the desired kernel

start = tic;
D = chunkermat(chnkr,fkern, opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)


A = [-1/2, 0; 1, -1/2];                                     % jump matrix (for exterior problem)

M = kron(eye(chnkr.npt), A);  

src = randn(2,1)*0.1;

[y1, grad, ~, ~, ~] = flex2d.hkdiffgreen(zk, src, chnkr.r);


nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);

normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC                          

firstbc = 1/(2*zk^2).*y1 ;
secondbc = 1/(2*zk^2).*normalderiv ;

[nt, ~] = size(D);
lhs = M + D;
%%
%

rhs = zeros(nt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;



tic
%sol = gmres(lhs, rhs, [], 1e-13, 200);
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = (-2:0.1:2)+0.01*randn();                                    % generate some targets
ys = (-2:0.1:2)+0.01*randn();
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

tic
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


[val, ~] = flex2d.hkdiffgreen(zk, src, targets(:,out));

trueval = 1/(2*zk^2).*val;

utarg(out) = Dsol;
true_sol(out) = trueval;


uerr = abs(utarg - true_sol)/(max(abs(true_sol(:))));
uerr = reshape(uerr,size(X));
figure
h = pcolor(X,Y,log10(uerr));
set(h,'EdgeColor','None'); hold on;
title("Absolute error (clamped plate kernel on a circle)", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

