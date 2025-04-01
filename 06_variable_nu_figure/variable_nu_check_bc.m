
% CREATING GEOMETRY 

zk = 4;
nus = [0.3; 0; -0.3];
utots = {};

cparams = [];
cparams.maxchunklen = 2/zk;       % setting a chunk length helps when the
                                    % frequency is known'

chnkr = chunkerfuncuni(@(t) cavity(t,10.5),64);
%chnkr = chunkerfunc(@(t) ellipse(t,3,1), cparams);
chnkr = chnkr.move([0;0],[0;0],0,2);
chnkr = sort(chnkr);
chnkr.npt

coefs = [nu; 0];
opts = [];

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
hold on
axis equal
drawnow

bdy_r = chnkr.r(:,8,10);
bdy_n = chnkr.n(:,8,10);
bdy_tau = chnkr.d(:,8,10);

h = 10^(-3);

d1f = [-49/20	6	-15/2	20/3	-15/4	6/5	-1/6 0 0] .';
d2f = [469/90	-223/10	879/20	-949/18	41	-201/10	1019/180	-7/10 0].';
d3f = [-801/80	349/6	-18353/120	2391/10	-1457/6	4891/30	-561/8	527/30	-469/240].';

d1f = [-3/2	2	-1/2 0 0].';
d2f = [2	-5	4	-1 0].';
d3f = [-5/2	9	-12	7	-3/2].';

% d1c = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280].';
% d2c = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560].';
% d3c = [-7/240	3/10	-169/120	61/30	0	-61/30	169/120	-3/10	7/240].';
% d4c = [7/240	-2/5	169/60	-122/15	91/8	-122/15	169/60	-2/5	7/240].';

d1c = [1/12	-2/3	0	2/3	-1/12].';
d2c = [-1/12	4/3 -5/2	4/3	-1/12].';
d3c = [-1/2	1	0	-1	1/2].';
d4c = [1	-4	6	-4	1].';

bilap = zeros(5);
bilap(:,3) = d4c;
bilap(3,:) = bilap(3,:) + d4c.';
bilap = bilap + 2*(d2c*d2c.');
bilap = bilap / h^4;

d3dndt2 = d1c*d2c.' / h^3;

d3dn3 = zeros(5);
d3dn3(:,3) = d3c / h^3;

d2dn2 = zeros(5);
d2dn2(:,3) = d2c / h^2;

d2dt2 = zeros(5);
d2dt2(3,:) = d2c / h^2;

kappa = signed_curvature(chnkr);
bdy_kappa = kappa(8,10);

theta = atan2(bdy_n(2),bdy_n(1)) - pi/2;

[tx1,tx2] = meshgrid(-2*h:h:2*h, 0:h:4*h);

R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Rotate and translate grid points
txy_rotated = R * [tx1(:)'; tx2(:)'];
txy_rotated(1,:) = txy_rotated(1,:) + bdy_r(1) + 1*h*bdy_n(1);
txy_rotated(2,:) = txy_rotated(2,:) + bdy_r(2) + 1*h*bdy_n(2);
targets = txy_rotated;

%tx_rot = txy_rotated(1, :) + bdy_r(1); %+ (4*h+1e-3)*bdy_n(1);
%ty_rot = txy_rotated(2, :) + bdy_r(2); %+ (4*h+1e-3)*bdy_n(2);

% Plot the transformed grid
plot(targets(1,:), targets(2,:), 'ko', 'MarkerFaceColor', 'k');
axis equal;
xlabel('X'); ylabel('Y');
grid on;
title('Rotated and Translated 9x9 Grid');
drawnow

%tx1 = reshape(tx_rot(:),size(tx1));
%tx2 = reshape(ty_rot(:),size(tx1));

%targets = [tx1(:) tx2(:)].';

targinfo = [];
targinfo.r = txy_rotated;
targinfo.n = repmat(chnkr.n(:,8,10),1,25) ;
targinfo.d = repmat(chnkr.d(:,8,10),1,25) ;
targinfo.d2 = repmat(chnkr.d2(:,8,10),1,25) ;

srcinfo = [];
srcinfo.r = chnkr.r(:,2,2);
srcinfo.n = chnkr.n(:,2,2);
srcinfo.d = chnkr.d(:,2,2);

fker1 =  flex2d.kern(zk, srcinfo, targinfo, 'free plate first part', coefs);        % build the desired kernel
eker1 =  flex2d.kern(zk, srcinfo, targinfo, 'free plate eval second', coefs);        % build the desired kernel

eker1 = reshape(eker1(:),[5 5]);

t1 = sum(d2dn2.*eker1,'all');
t2 = nu*sum(d2dt2.*eker1,'all');
t3 = sum(d3dn3.*eker1,'all');
t4 = (2-nu)*sum(d3dndt2.*eker1,'all');
t5 = (1-nu)*bdy_kappa*(sum(d2dt2.*eker1,'all'));
t6 = (1-nu)*bdy_kappa*(-sum(d2dn2.*eker1,'all'));

bc1 = t1 + t2;
bc2 = t3 + t4 + t5 + t6;

targinfo2 = [];
targinfo2.r = bdy_r;
targinfo2.n = chnkr.n(:,8,10);
targinfo2.d = chnkr.d(:,8,10);
targinfo2.d2 = chnkr.d2(:,8,10);
fker2 =  flex2d.kern(zk, srcinfo, targinfo2, 'free plate first part', coefs);        % build the desired kernel

err1 = abs(- fker2(1,2) + bc1) / max(abs(t1),abs(t2))
err2 = abs(- fker2(2,2) + bc2) / max([abs(t3),abs(t4),abs(t5),abs(t6)])
% err3 = abs(sum(bilap.*eker1,'all')-zk^4*eker1(3,3)) 

% sum(d2dn2.*eker1,'all')+nu*sum(d2dt2.*eker1,'all') - fker1(25,2)
% 
% sum(d3dn3.*eker1,'all')+(2-nu)*sum(d3dndt2.*eker1,'all')+...
%     +(1-nu)*bdy_kappa*(sum(d2dt2.*eker1-d2dn2.*eker1,'all')) - fker1(26,2)


% d3dndt2 = d1f*d2c.' / h^3;
% 
% d3dn3 = zeros(5);
% d3dn3(:,3) = d3f / h^3;
% 
% d2dn2 = zeros(5);
% d2dn2(:,3) = d2f / h^2;
% 
% d2dt2 = zeros(5);
% d2dt2(1,:) = d2c / h^2;
% 
% t1 = sum(d2dn2.*eker1,'all');
% t2 = nu*sum(d2dt2.*eker1,'all');
% t3 = sum(d3dn3.*eker1,'all');
% t4 = (2-nu)*sum(d3dndt2.*eker1,'all');
% t5 = (1-nu)*bdy_kappa*(sum(d2dt2.*eker1,'all'));
% t6 = (1-nu)*bdy_kappa*(-sum(d2dn2.*eker1,'all'));
% 
% bc1 = t1 + t2;
% bc2 = t3 + t4 + t5 + t6;
% 
% err1 = abs(- fker2(1,2) + bc1) / max(abs(t1),abs(t2))
% err2 = abs(- fker2(2,2) + bc2) / max(abs([t3,t4,t5,t6]))


for ii = 1:3
    
    nu = nus(ii);

    kappa = signed_curvature(chnkr);
    kappa = kappa(:)';

    theta = -pi;
    d = -[cos(theta) sin(theta)];
    
    [val, grad, hess, third] = planewave1(zk, chnkr.r(:,:), d);
    
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
    
    % thetas = -pi:0.05:pi;
    % R=1;
    % xs = R*cos(thetas);
    % ys = R*sin(thetas);
    % targets = [xs; ys];
    % [~,na] = size(targets);
    
    ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
    ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
    ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);
    
    coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets);
    
    start1 = tic;
    Dsol = chunkerkerneval(chnkr, ikern1,rho1, targets) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets);
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)
    
    Dsol = reshape(Dsol,size(tx1));

    t1 = sum(d2dn2.*Dsol,'all');
    t2 = nu*sum(d2dt2.*Dsol,'all');
    t3 = sum(d3dn3.*Dsol,'all');
    t4 = (2-nu)*sum(d3dndt2.*Dsol,'all');
    t5 = (1-nu)*bdy_kappa*(sum(d2dt2.*Dsol,'all'));
    t6 = (1-nu)*bdy_kappa*(-sum(d2dn2.*Dsol,'all'));

    bc1 = t1 + t2;
    bc2 = t3 + t4 + t5 + t6;

    err1 = abs(firstbc(152) + bc1) / max(abs(t1),abs(t2))
    err2 = abs(secondbc(152) + bc2) / max([abs(t3),abs(t4),abs(t5),abs(t6)])

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
% clim([minval maxval])
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
% clim([minval maxval])
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

minval = max([min(real(utots{1}(:))) min(real(utots{2}(:))) min(real(utots{3}(:)))]);
maxval = min([max(real(utots{1}(:))) max(real(utots{2}(:))) max(real(utots{3}(:)))]);


figure(3)
tiledlayout(1,3,"TileSpacing","compact")
nexttile
h = pcolor(Tx,Ty,real(utots{1}));
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
h = pcolor(Tx,Ty,real(utots{2}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
% clim([minval maxval])
plot(chnkr,'k-','LineWidth',2);
title(['\nu = ',num2str(nus(2))])
axis square

ax = gca;  % Get the current axis
ax.XAxis.LineWidth = 0.8;  % Set the X-axis tick mark width
ax.YAxis.LineWidth = 0.8;  % Set the Y-axis tick mark width
set(ax, 'FontSize',12)

nexttile
h = pcolor(Tx,Ty,real(utots{3}));
h.EdgeColor = 'none';
h.FaceColor = 'interp';
hold on
colorbar
% clim([minval maxval])
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