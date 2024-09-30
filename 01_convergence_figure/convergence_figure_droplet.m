clear 
addpath(genpath(pwd)) 

maxchunklens = [ 2.7 2.5  1.8  1.2  0.8 0.75 0.7 ]; 
npts = maxchunklens*0;
free_errors = maxchunklens*0;
clamped_errors = maxchunklens*0;
supported_errors = maxchunklens*0;

zk = 8;
nu = 1/3;
cparams = [];

cparams.eps = 1e-5;
thetas = 0:pi/6:2*pi-pi/6;
[targets, ~, ~] = droplet(thetas);
targets = targets*1.5;
%targets = [3*cos(thetas); 1.5*sin(thetas)];
centre = [0.8 ; 0.5];

for i = 1:numel(npts)

    
    length = maxchunklens(i);
    cparams.maxchunklen = length;       % setting a chunk length helps when the
                                        % frequency is known'
    chnkr = chunkerfunc(@(t) droplet(t), cparams);
    chnkr = chnkr.sort();
    npts(i) = chnkr.npt;

    figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
    clf
    plot(chnkr, '-x')
    hold on
    quiver(chnkr)
    hold on 
    scatter(targets(1,:),targets(2,:),36,'filled')
    hold on 
    scatter(centre(1),centre(2),36,'filled')
    axis equal
    drawnow
        
    coefs = [nu; 0];
    opts = [];
    opts.sing = 'log';

    kappa = signed_curvature(chnkr);
    kappa = kappa(:);

    % clamped plate
    
    opts = [];
    opts.sing = 'log';

    opts2 = [];
    opts2.quad = 'native';
    opts2.sing = 'smooth';

    fkern =  @(s,t) flex2d.kern(zk, s, t, 'clamped-plate');           % build the desired kernel
    fkern2 =  @(s,t) flex2d.kern(zk, s, t, 'clamped-plate K21');           % build the desired kernel

    start = tic;
    D = chunkermat(chnkr,fkern, opts);
    M2 = chunkermat(chnkr,fkern2, opts2);
    M2(isnan(M2)) = 0;

    M2 = M2 - diag(3*kappa.^2/(4*pi).*chnkr.wts(:));
    M2 = 0;

    D(1:2:end,1:2:end) = D(1:2:end,1:2:end) - 0.5*eye(chnkr.npt);
    D(2:2:end,1:2:end) = D(2:2:end,1:2:end) + M2 + kappa.*eye(chnkr.npt);
    D(2:2:end,2:2:end) = D(2:2:end,2:2:end) - 0.5*eye(chnkr.npt);


    t1 = toc(start);
    fprintf('%5.2e s : time to assemble matrix\n',t1)

    [y1, grad, ~, ~, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);

    nx = chnkr.n(1,:); 
    ny = chnkr.n(2,:);

    normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         

    firstbc = 1/(2*zk^2).*y1;
    secondbc = 1/(2*zk^2).*normalderiv;


    rhs = zeros(2*chnkr.npt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;


    sol = D\rhs;

    rho1 = sol(1:2:end);                                    % first density
    rho2 = sol(2:2:end);        


    ikern1 = @(s,t) flex2d.kern(zk, s, t, 'first kernel');                              % build the kernel of evaluation          
    ikern2 = @(s,t) flex2d.kern(zk, s, t, 'second kernel');


    start1 = tic;
    utarg = chunkerkerneval(chnkr, ikern1,rho1, targets) + chunkerkerneval(chnkr, ikern2, rho2, targets);
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


    [val, ~] = flex2d.hkdiffgreen(zk, centre, targets);

    true_sol = 1/(2*zk^2).*val;

    uerr = utarg - true_sol;
%    uerr = uerr ./  max(abs(true_sol));
    uerr = uerr ./  (chnkr.wts(:).'*(abs(rho1) + abs(rho2)));
    clamped_errors(i) = max(abs(uerr));

    % free plate

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

    opts3 = [];
    opts3.sing = 'pv';

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

    zkimag = (1i)*zk;
    [~, ~, hess, third, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);

    nx = chnkr.n(1,:).'; 
    ny = chnkr.n(2,:).';

    dx = chnkr.d(1,:).';
    dy = chnkr.d(2,:).';

    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds); % normalization
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

    [nt, ~] = size(sysmat);

    rhs = zeros(nt, 1); 
    rhs(1:2:end) = firstbc ; 
    rhs(2:2:end) = secondbc;

    sol = lhs\rhs;

    rho1 = sol(1:2:end);                                    % first density
    rho2 = sol(2:2:end);        


    ikern1 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
    ikern2 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval second');
    ikern3 = @(s,t) flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

    coupled = chunkerkerneval(chnkr, ikern3, H*rho1, targets);


    start1 = tic;
    utarg = chunkerkerneval(chnkr, ikern1,rho1, targets) + coupled +chunkerkerneval(chnkr, ikern2, rho2, targets);
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)



    [val,~] = flex2d.hkdiffgreen(zk,centre,targets);        % Hankel part

    true_sol = 1/(2*zk^2).*val ;

    uerr = utarg - true_sol;
%    uerr = uerr ./  max(abs(true_sol));
    uerr = uerr ./  (chnkr.wts(:).'*(abs(rho1) + abs(rho2)));
    free_errors(i) = max(abs(uerr));

    % supported plate

    ikern1 =  @(s,t) flex2d.suppkern(zk, s, t, 'supported plate droplet',coefs);           % build the desired kernel
    ikern2 =  @(s,t) flex2d.suppkern(zk, s, t, 'supported plate K21 droplet',coefs);           % build the desired kernel

    opts = [];
    opts.sing = 'log';
    
    opts2 = [];
    opts2.quad = 'native';
    opts2.sing = 'smooth';
    
    start = tic;
    M = chunkermat(chnkr,ikern1, opts);
    M2 = chunkermat(chnkr,ikern2, opts2);
    M2(isnan(M2)) = 0;
    
    % droplet: 

    xs = chnkr.r(1,:,:);
    ys = chnkr.r(2,:,:);

    xs = xs(:);
    ys = ys(:);

    a = 2;
    b = 1;
    c = -0.4;

    t = atan2(a*ys - c*xs.^2/a, b*xs);

    % x = a*cos(t);
    % y = b*sin(t) + c*cos(t).^2;
    
    x1 = -a*sin(t);
    y1 = b*cos(t) - 2*c*cos(t).*sin(t);
    
    x2 = -a*cos(t);
    y2 = -b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;

    x3 = a*sin(t);
    y3 = -b*cos(t) + 8*c*sin(t).*cos(t);

    x4 = a*cos(t);
    y4 = b*sin(t) + 8*c*cos(t).^2 - 8*c*sin(t).^2;

    kp2 = (-x4.*y1.^3 + y1.^2.*(4*y3.*x2+y4.*x1) + x1.*(-3*x2.^2.*y2 + y4.*x1.^2 - 4*x3.*x1.*y2 - 3*y2.^3) ...
         + y1.*(3*x2.^3 - x1.*(x4.*x1 + 4*y3.*y2) + x2.*(4*x3.*x1 + 3*y2.^2))) ./ (x1.^2 + y1.^2).^(7/2) ...
         - 6*(x1.*x2 + y1.*y2).*((x1.^2+y1.^2).*(x1.*y3 - y1.*x3) + 3*(y1.*x2 - x1.*y2).*(x1.*x2 + y1.*y2)) ./ (x1.^2 + y1.^2).^(9/2);
    
    k21diag = (nu - 1)*(12*kappa.^3*(nu^2 - nu + 4) + kp2*(-5*nu^2 + 4*nu + 33))/(48*pi*(nu - 3)) ; %+ 1i*zk^2/64*(3+nu)*((-1+nu)*(7+nu)/(3-nu))*kappa;
    
    K21fo = M2 + diag(k21diag.*chnkr.wts(:));
    
    c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

    M(1:2:end,1:2:end) = M(1:2:end,1:2:end) - 0.5*eye(chnkr.npt) ;
    M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + K21fo + c0.*kappa.^2.*eye(chnkr.npt);
    M(2:2:end,2:2:end) = M(2:2:end,2:2:end) - 0.5*eye(chnkr.npt);
    
    
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
    
    [nt, ~] = size(M);
    
    rhs = zeros(nt, 1); 
    rhs(1:2:end) = firstbc ; 
    rhs(2:2:end) = secondbc;
    
    start = tic;
    sol = M\rhs;
    t3 = toc(start); 
    fprintf('%5.2e s : time to solve system \n',t3);

    
    rho1 = sol(1:2:end);                                    % first density
    rho2 = sol(2:2:end);  % second density
    
    
    ikern1 = @(s,t) flex2d.suppkern(zk, s, t, 'supported plate K1 eval droplet',coefs);                              % build the kernel of evaluation          
    ikern2 = @(s,t) flex2d.suppkern(zk, s, t, 'supported plate K2 eval',coefs);
    
    start1 = tic;
    utarg = chunkerkerneval(chnkr, ikern1, rho1, targets) + ...
            chunkerkerneval(chnkr, ikern2, rho2, targets);

    [val, ~] = flex2d.hkdiffgreen(zk, centre, targets);
    
    true_sol = 1/(2*zk^2).*val;

    uerr = utarg - true_sol;
%    uerr = uerr ./  max(abs(true_sol));
    uerr = uerr ./  (chnkr.wts(:).'*(abs(rho1) + abs(rho2)));
    supported_errors(i) = max(abs(uerr));
end

%%

figure(2)
tiledlayout(1,2)
nexttile
loglog(npts , clamped_errors, '.-', npts, free_errors,'.-', npts, supported_errors,'.-', MarkerSize=20)
%loglog(npts, supported_errors,'.-', MarkerSize=20)
hold on 
loglog(npts(2:end-2), 10^26.5.*npts(2:end-2).^(-16),'--k')
hold on
%legend('supported plate','n^{-12}')
legend('Clamped Plate','Free Plate','Supported Plate', 'n^{-16}')
xlabel('N')
ylabel('Relative error')
title('Analytic solution test')

nexttile
length = maxchunklens(1);
cparams.maxchunklen = length;       % setting a chunk length helps when the
                                    % frequency is known'
chnkr = chunkerfunc(@(t) droplet(t), cparams);
chnkr = chnkr.sort();

plot(chnkr, '-xk')
hold on 
scatter(targets(1,:),targets(2,:),36,'filled')
hold on 
scatter(centre(1),centre(2),36,'filled')
xlim([-3.5 3.5])
ylim([-2.5 2.5])

% Integral of absolute value of density on boundary 
% Maximum of density on boundary 


rmpath(genpath(pwd)) 
