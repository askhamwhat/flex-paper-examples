
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


xs = -1:0.01:2;                                    % generate some targets
ys = -2:0.01:1;
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
set(gca,'FontName','CMU Serif')
set(gca,'FontSize',20)
clim([-1,1]);
colorbar

nexttile
h = pcolor(X,Y,real(us));
h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
title("Scattered field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
set(gca,'FontName','CMU Serif')
set(gca,'FontSize',20)
clim([-1,1]);
colorbar

nexttile
h = pcolor(X,Y,real(utot));
h.FaceColor = 'interp';
set(h,'EdgeColor','None'); hold on;
title("Total field", 'FontSize',16)
plot(chnkr,'w-','LineWidth',2);
set(gca,'FontName','CMU Serif')
set(gca,'FontSize',20)
clim([-1,1]);
colorbar


