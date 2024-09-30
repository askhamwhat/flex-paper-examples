src = [2;-3];
%src = [2;0];
[~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, src, chnkr.r);
[~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), src, chnkr.r);

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


xs = 0:0.1:4;                                    % generate some targets
ys = -6:0.1:0;
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

[val,~] = chnk.flex2d.green(zk,src,targets(:,out));        % Hankel part

zkimag = (1i)*zk;
[valK,~] = chnk.flex2d.green(zkimag,src, targets(:,out));    % modified bessel K part

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
title("Computed Field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile
h = pcolor(X,Y, real(true_sol));
set(h,'EdgeColor','None'); hold on;
title("True field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar

nexttile 
h = pcolor(X,Y,log10(abs(uerr)./max(abs(true_sol(:)))));
set(h,'EdgeColor','None'); hold on;
title("Log 10 relative error", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
colorbar
