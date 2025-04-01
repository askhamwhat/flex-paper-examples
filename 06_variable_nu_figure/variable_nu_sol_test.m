
% CREATING GEOMETRY 

zk = 12;
nus = [0.3; 0; -0.3];
utots = {};

cparams = [];
cparams.maxchunklen = 4/zk;       % setting a chunk length helps when the
                                    % frequency is known'

chnkr = chunkerfuncuni(@(t) cavity(t,10.5),64);
%chnkr = chunkerfunc(@(t) ellipse(t,3,1), cparams);
chnkr = chnkr.move([0;0],[0;0],0,2);
chnkr.npt

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal
drawnow

centre = [1.1; 0];

[~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds); % normalization
tauy = (dy./ds);

kappa = signed_curvature(chnkr);
kappa = kappa(:);

h = 0.05;
[Tx,Ty] = meshgrid(-3:h:3);
targets = [Tx(:) Ty(:)].';

in = chunkerinterior(chnkr, targets); 
out = ~in; 

[na,~] = size(targets);
utarg = zeros(na, 1);


for ii = 1:3
    
    nu = nus(ii);
    
    coefs = [nu; 0];
    opts = [];
    opts.sing = 'log';
    
    opts2 = [];
    opts2.quad = 'native';
    opts2.sing = 'smooth';
    
    opts3 = [];
    opts3.sing = 'pv';
    
    
    fkern1 =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel
    fkern1bh =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part bh', coefs);        % build the desired kernel
    fkern2bh =  @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert bh', coefs);        % build the desired kernel
    fkern2 =  @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert', coefs);        % build the desired kernel
    double = @(s,t) lap2d.kern(s,t,'d',coefs);
    hilbert = @(s,t) lap2d.kern(s,t,'hilb',coefs);
    
    sysmat1 = chunkermat(chnkr,fkern1, opts);
    sysmat1bh = chunkermat(chnkr,fkern1bh, opts2);
    sysmat2bh = chunkermat(chnkr,fkern2bh, opts2);
    sysmat2 = chunkermat(chnkr,fkern2, opts);
    
    D = chunkermat(chnkr, double, opts);
    
    H = chunkermat(chnkr, hilbert, opts3);     
    
    % Perform diagonal replacement for smooth quads here
    
    sysmat1bh(isnan(sysmat1bh)) = 0;
    sysmat1bh(2:2:end,1:2:end) = sysmat1bh(2:2:end,1:2:end) + diag((-3+3*nu)/(8*pi)*kappa.^2.*chnkr.wts(:));
    sysmat1 = sysmat1 + sysmat1bh;
    
    sysmat2bh(isnan(sysmat2bh)) = 0;
    sysmat2 = sysmat2 + sysmat2bh;
    
    sysmat2(1:2:end,1:2:end) = sysmat2(1:2:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
    sysmat2(2:2:end,1:2:end) = sysmat2(2:2:end,1:2:end)*H;
    
    sysmat = sysmat1 + sysmat2 ;
    
    D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
    D = kron(eye(chnkr.npt), D);
    
    lhs =  D + sysmat;
    % cond(lhs)
    
    nx = chnkr.n(1,:).'; 
    ny = chnkr.n(2,:).';
    
    dx = chnkr.d(1,:).';
    dy = chnkr.d(2,:).';
    
    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds);                                                                       % normalization
    tauy = (dy./ds);

    firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny))+...
           coefs(1)/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy));


    secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
           third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))+...
            (2-coefs(1))/(2*zk^2).*(third(:, :, 1).*(taux.*taux.*nx) + third(:, :, 2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
            third(:, :, 3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
            + third(:, :, 4).*(tauy.*tauy.*ny))+...
            (1-coefs(1)).*(kappa).*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
            (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)));


    rhs = zeros(2*chnkr.npt, 1); 
    rhs(1:2:end) = firstbc ; 
    rhs(2:2:end) = secondbc;
   
    
    tic
    sol = lhs\rhs;
    toc;
    
    
    rho1 = sol(1:2:end);                                    % first density
    rho2 = sol(2:2:end);        
    
    % thetas = -pi:0.05:pi;
    % R=1;
    % xs = R*cos(thetas);
    % ys = R*sin(thetas);
    % targets = [xs; ys];
    % [~,na] = size(targets);
    
    ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
    ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
    ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);
    
    coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets(:, out));
    
    start1 = tic;
    utarg(out) = chunkerkerneval(chnkr, ikern1,rho1, targets(:, out)) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets(:,out));
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)
    
    true_sol = 0*utarg;
    true_sol(out) = flex2d.hkdiffgreen(zk,centre,targets(:,out));        % Hankel part

    true_sol = 1/(2*zk^2).*true_sol ;

    uerr = utarg - true_sol;
%    uerr = uerr ./  max(abs(true_sol));
    uerr = uerr ./  (chnkr.wts(:).'*(abs(rho1) + abs(rho2)));
    uerr = reshape(uerr,size(Tx));
    uerr(uerr == 0) = NaN;
    utots{ii} = abs(uerr);
    usols{ii} = reshape(utarg,size(Tx));

end

%%

% figure(1)
% tiledlayout(1,2,"TileSpacing","compact")
% nexttile
% h = pcolor(Tx,Ty,abs((utots{1})-(utots{2})));
% h.EdgeColor = 'none';
% h.FaceColor = 'interp';
% hold on
% colorbar
% %clim([-max(abs(us{1}(:)-us{2}(:))) max(abs(us{1}(:)-us{2}(:)))])
% plot(chnkr,'k-','LineWidth',2);
% axis square 
% title(['\delta\nu = ',num2str(nus(1) - nus(2)),' (\nu_{ice} - \nu_{cork})'])
% 
% ax = gca;  % Get the current axis
% ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
% ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
% set(ax, 'FontSize',12)
% 
% nexttile
% h = pcolor(Tx,Ty,(abs(utots{3}-utots{2})));
% h.EdgeColor = 'none';
% h.FaceColor = 'interp';
% hold on
% colorbar
% %clim([-max(abs(us{3}(:)-us{2}(:))) max(abs(us{3}(:)-us{2}(:)))])
% plot(chnkr,'k-','LineWidth',2);
% title(['\delta\nu = ',num2str(nus(3) - nus(2)), ' (\nu_{auxetic} - \nu_{cork})'])
% axis square
% 
% ax = gca;  % Get the current axis
% ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
% ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
% set(ax, 'FontSize',12)
% 
% 
% fontname(gcf, 'Helvetica')
% 
% set(gcf,'Position',[541 592 700 317])

%%
% saveas(gcf,'variable_nu.fig','fig')
% exportgraphics(gcf,'variable_nu.pdf','ContentType','vector')


%% plotting absolute value of the total field

minval = min([min(abs(utots{1}(:))) min(abs(utots{2}(:))) min(abs(utots{3}(:)))]);
maxval = max([max(abs(utots{1}(:))) max(abs(utots{2}(:))) max(abs(utots{3}(:)))]);


figure(2)
tiledlayout(1,3,"TileSpacing","compact")
nexttile
h = pcolor(Tx,Ty,abs(utots{1}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
% clim([minval maxval])
plot(chnkr,'k-','LineWidth',2);
axis square 
title(['\nu = ',num2str(nus(1))])

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

nexttile
h = pcolor(Tx,Ty,abs(utots{2}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
plot(chnkr,'k-','LineWidth',2);
title(['\nu = ',num2str(nus(2))])
axis square

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

nexttile
h = pcolor(Tx,Ty,abs(utots{3}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
plot(chnkr,'k-','LineWidth',2);
title(['\nu = ',num2str(nus(3))])
axis square

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

fontname(gcf, 'Helvetica')

set(gcf,'Position',[541 592 927 317])


% plotting real part of solution instead of absolute value



figure(3)
tiledlayout(1,3,"TileSpacing","compact")
nexttile
h = pcolor(Tx,Ty,real(usols{1}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
plot(chnkr,'k-','LineWidth',2);
axis square 
title(['\nu = ',num2str(nus(1))])

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

nexttile
h = pcolor(Tx,Ty,real(usols{2}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
plot(chnkr,'k-','LineWidth',2);
title(['\nu = ',num2str(nus(2))])
axis square

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

nexttile
h = pcolor(Tx,Ty,real(usols{3}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
plot(chnkr,'k-','LineWidth',2);
title(['\nu = ',num2str(nus(3))])
axis square

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

fontname(gcf, 'Helvetica')

set(gcf,'Position',[541 592 927 317])

return