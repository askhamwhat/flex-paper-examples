%clear
%clc

zk = 10;
nu = 1/3;

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
pref = [];
pref.k = 16;
cparams.maxchunklen = 2*pi / zk; % setting a chunk length helps when the
                              % frequency is known
cparams.splitatpoints = true;
cparams.pref = pref;
pos = readmatrix('penguin_pos_2.txt');

% rlft = pos(1:(end-2),:);
% rrgt = pos(3:(end)  ,:);
% rcen = pos(2:(end-1),:);
% d1 = rlft-rcen;
% n1 = vecnorm(d1,2,2);
% d1 = d1./n1;
% d2 = rrgt-rcen;
% n2 = vecnorm(d2,2,2);
% d2 = d2./n2;
% csgn = sum(d1.*d2,2);
% 
% inds =find(csgn >-0.95);
% numel(inds)
% inds = [1;inds+1;size(pos,1)];
% inds = sort(inds);
%pos = pos(inds,:)/100;
pos = pos/100;
chnkr = chunkerfit(pos.',cparams);

%chnkr = chunkerpoly(pos.', cparams);
%chnkr = refine(chnkr,cparams);

% 
% pref = []; 
% pref.k = 16;
% narms = 5;
% amp = 0.25;
% start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
% t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
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

x = chnkr.d(1,:);
y = chnkr.d(2,:);

denom = sqrt(x.^2+y.^2);
numer = chnkr.d(1,:).*chnkr.d2(2,:)-chnkr.d2(1,:).*chnkr.d(2,:);

kappa = numer./(denom.^3);
kappa = signed_curvature(chnkr(:));
kappa = kappa(:).';
size(kappa)
hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
hilb2 = sysmat2*H ;

mat1 =  sysmat(1:2:end, 1:2:end);
mat4 = sysmat(2:2:end, 2:2:end);



sysmat(1:2:end, 1:2:end) = mat1 + hilb;
sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa.').*(K21hilbert*H + K21second);
sysmat(2:2:end, 2:2:end) = mat4 + (kappa.').*(K22second);


A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                    

M = kron(eye(chnkr.npt), A);

lhs =  M + sysmat ;

[~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, [0;0], chnkr.r);
[~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), [0;0], chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

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
        (1-coefs(1)).*(kappa.').*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
         1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
        hessK(:, :, 3).*tauy.*tauy)-...
        (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
       1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));



[nt, ~] = size(sysmat);

rhs = zeros(nt, 1); 
rhs(1:2:end) = firstbc ; 
rhs(2:2:end) = secondbc;

tic
sol = lhs\rhs;
toc;

rho1 = sol(1:2:end);                                    % first density
rho2 = sol(2:2:end);        


xs = -2:0.1:2;                                    % generate some targets
ys = -2:0.1:2;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];


% rmin = min(chnkr); rmax = max(chnkr);                                       % plot commands
% xl = rmax(1)-rmin(1);
% yl = rmax(2)-rmin(2);
% nplot = 200;
% xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot); 
% ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
% [xxtarg,yytarg] = meshgrid(xtarg,ytarg);
% targets = zeros(2,length(xxtarg(:)));
% targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);
[~,na] = size(targets);



tic
in = chunkerinterior(chnkr, targets); 
out = ~in; 
toc


ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));

start1 = tic;
Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

[val,~] = chnk.flex2d.green(zk,[0;0],targets(:,out));        % Hankel part

zkimag = (1i)*zk;
[valK,~] = chnk.flex2d.green(zkimag,[0;0], targets(:,out));    % modified bessel K part

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;


utarg(out) = Dsol;
true_sol(out) = trueval;

utarg = reshape(utarg,size(X));
true_sol = reshape(true_sol,size(X));
uerr = utarg - true_sol;
figure
tiledlayout(1,3)
nexttile 
h = pcolor(X,Y, real(utarg));
set(h,'EdgeColor','None'); hold on;
title("BIE Solution", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile
h = pcolor(X,Y, real(true_sol));
set(h,'EdgeColor','None'); hold on;
title("True Solution", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile 
h = pcolor(X,Y,log10(abs(uerr) / max(abs(true_sol(:)))));
set(h,'EdgeColor','None'); hold on;
title("Relative error (free plate kernel)", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

