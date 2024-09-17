clear 


% zk = 1;
% 
% %
% 
% zk = 0.1;               % our k (wave number)

maxchunklens = [4 2 1 0.5 0.25 0.125];
npts = maxchunklens*0;
free_errors = maxchunklens*0;
clamped_errors = maxchunklens*0;
supported_errors = maxchunklens*0;

zk = 5;
nu = 1/3;
cparams = [];

cparams.eps = 1e-5;
cparams.nover = 0;
thetas = 0:pi/12:2*pi-pi/12;
targets = [3*cos(thetas); 1.5*sin(thetas)];
centre = [0.5;0.5];

for i = 1:numel(npts)

    length = maxchunklens(i);
    cparams.maxchunklen = length;       % setting a chunk length helps when the
                                        % frequency is known'
    chnkr = chunkerfunc(@(t) ellipse(t), cparams);
    npts(i) = chnkr.npt;

    figure(1)                                                   % plot the chunker-object (supposed to be a circle centered at 1 with radius 1)
    clf
    plot(chnkr, '-x')
    hold on
    quiver(chnkr)
    hold on 
    %scatter(xs(:),ys(:),36,kp,'filled')
    axis equal
    drawnow
        
    coefs = [nu; 0];
    opts = [];
    opts.sing = 'log';

    % clamped plate
    
    opts = [];
    opts.sing = 'log';
    fkern =  @(s,t) flex2d.kern(zk, s, t, 'clamped-plate');           % build the desired kernel
    
    start = tic;
    D = chunkermat(chnkr,fkern, opts);
    t1 = toc(start);
    
    fprintf('%5.2e s : time to assemble matrix\n',t1)
    
    
    
    kappa = signed_curvature(chnkr);
    kappa = kappa(:);
    
    A = zeros(2, 2, chnkr.npt);
    start = tic;
    for j = 1:chnkr.npt
        A(:, :, j) = [-0.5, 0 ; kappa(j), -0.5];
    end
    t3 = toc(start); 
    fprintf('%5.2e s : time to construct jump matrix\n',t3);
    
    K = num2cell(A, [1 2]);
    M = blkdiag(K{:}); 
     
    
    [y1, grad, ~, ~, ~] = flex2d.hkdiffgreen(zk, centre, chnkr.r);
    
    nx = chnkr.n(1,:); 
    ny = chnkr.n(2,:);
    
    normalderiv = grad(:, :, 1).*(nx.')+ grad(:, :, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         
    
    firstbc = 1/(2*zk^2).*y1;
    secondbc = 1/(2*zk^2).*normalderiv;
    
    [nt, ~] = size(D);
    lhs = M + D;
    
    rhs = zeros(nt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;
    
    
    tic
    %sol = gmres(lhs, rhs, [], 1e-13, 400);
    sol = lhs\rhs;
    toc;
    
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
    uerr = uerr ./  abs(true_sol);
    clamped_errors(i) = norm(uerr);

    % free plate

    fkern =  @(s,t) flex2d.kern(zk, s, t, 'free plate first part', coefs);        % build the desired kernel
    
    
    fkern1 = @(s,t) flex2d.kern(zk, s, t, 'free plate hilbert subtract', coefs);                   % hilbert subtraction kernels in K11
    fkern2 = @(s,t) flex2d.kern(zk, s, t, 'free plate coupled hilbert', coefs);   
    
    hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');
    double = @(s,t) chnk.lap2d.kern(s,t, 'd');
    
    fkern3 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 first part', coefs);                     % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)
    
    fkern4 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 second part', coefs);                    % kernels in K21 needs to multiply by curvature
    
    fkern5 = @(s,t) flex2d.kern(zk, s, t, 'free plate K21 hilbert part', coefs);                   % kernels in K21 coupled with hilbert transforms and needs to multiply by curvature
    
    fkern6 = @(s,t) flex2d.kern(zk, s, t, 'free plate K22 second part', coefs);                    % kernels in K22 needs to multiply by curvature
    
    
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
    
    
    
    kappa = signed_curvature(chnkr);
    kappa = kappa(:);
    
    
    hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
    hilb2 = sysmat2*H ;
    
    mat1 =  sysmat(1:2:end, 1:2:end);
    mat4 =  sysmat(2:2:end, 2:2:end);
    
    
    sysmat(1:2:end, 1:2:end) = mat1 + hilb;
    sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa).*(K21hilbert*H + K21second);
    sysmat(2:2:end, 2:2:end) = mat4 + (kappa).*(K22second);
    
    
    A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
    
    M = kron(eye(chnkr.npt), A);
    
    lhs =  M + sysmat;
    
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
    
    tic
    %sol = gmres(lhs, rhs, [], 1e-13, 200);
    sol = lhs\rhs;
    toc;
    
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
    uerr = uerr ./  abs(true_sol);
    free_errors(i) = norm(uerr);

    % supported plate

    ikern1 =  @(s,t) suppkern(zk, s, t, 'supported plate',coefs);           % build the desired kernel
    ikern2 =  @(s,t) suppkern(zk, s, t, 'supported plate K21',coefs);           % build the desired kernel

    opts = [];
    opts.sing = 'log';
    
    opts2 = [];
    opts2.quad = 'native';
    opts2.sing = 'smooth';
    
    start = tic;
    M = chunkermat(chnkr,ikern1, opts);
    M2 = chunkermat(chnkr,ikern2, opts2);
    M2(isnan(M2)) = 0;
    
    K11fo = M(1:2:end,1:2:end);
    K12fo = M(1:2:end,2:2:end);
    K21fo = M(2:2:end,1:2:end);
    K22fo = M(2:2:end,2:2:end);
    
    % Ellipse:  
    
    c = 2;
    
    xs = chnkr.r(1,:,:);
    ys = chnkr.r(2,:,:);
    
    xs = xs(:);
    ys = ys(:);
    
    t = atan2(c*ys,xs);
    
    x1 = -c*sin(t);    
    y1 = cos(t);
    
    x2 = -c*cos(t);
    y2 = -sin(t);
    
    x3 = c*sin(t);
    y3 = -cos(t);
    
    x4 = c*cos(t);
    y4 = sin(t);

    kp2 = (-x4.*y1.^3 + y1.^2.*(4*y3.*x2+y4.*x1) + x1.*(-3*x2.^2.*y2 + y4.*x1.^2 - 4*x3.*x1.*y2 - 3*y2.^3) ...
         + y1.*(3*x2.^3 - x1.*(x4.*x1 + 4*y3.*y2) + x2.*(4*x3.*x1 + 3*y2.^2))) ./ (x1.^2 + y1.^2).^(7/2) ...
         - 6*(x1.*x2 + y1.*y2).*((x1.^2+y1.^2).*(x1.*y3 - y1.*x3) + 3*(y1.*x2 - x1.*y2).*(x1.*x2 + y1.*y2)) ./ (x1.^2 + y1.^2).^(9/2);
    
    k21diag = (nu - 1)*(12*kappa.^3*(nu^2 - nu + 4) + kp2*(-5*nu^2 + 4*nu + 33))/(48*pi*(nu - 3)) + 1i*zk^2/64*(3+nu)*((-1+nu)*(7+nu)/(3-nu))*kappa;
    
    K21fo = M2 + diag(k21diag.*chnkr.wts(:));
    
    c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

    M(1:2:end,1:2:end) = K11fo - 0.5*eye(chnkr.npt) ;
    M(2:2:end,1:2:end) = K21fo + c0.*kappa.^2.*eye(chnkr.npt);
    M(1:2:end,2:2:end) = K12fo;
    M(2:2:end,2:2:end) = K22fo - 0.5*eye(chnkr.npt);
    
    lhs = M;
    
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
    
    [nt, ~] = size(lhs);
    
    rhs = zeros(nt, 1); 
    rhs(1:2:end) = firstbc ; 
    rhs(2:2:end) = secondbc;
    
    start = tic;
    sol = lhs\rhs;
    t3 = toc(start); 
    fprintf('%5.2e s : time to solve system \n',t3);

    
    rho1 = sol(1:2:end);                                    % first density
    rho2 = sol(2:2:end);  % second density
    
    
    ikern1 = @(s,t) suppkern(zk, s, t, 'supported plate K1 eval',coefs);                              % build the kernel of evaluation          
    ikern2 = @(s,t) suppkern(zk, s, t, 'supported plate K2 eval',coefs);
    
    start1 = tic;
    utarg = chunkerkerneval(chnkr, ikern1, rho1, targets) + ...
            chunkerkerneval(chnkr, ikern2, rho2, targets);

    [val, ~] = flex2d.hkdiffgreen(zk, centre, targets);
    
    true_sol = 1/(2*zk^2).*val;

    uerr = utarg - true_sol;
    uerr = uerr ./  abs(true_sol);
    supported_errors(i) = norm(uerr);
end

figure(2)
loglog(npts , clamped_errors, npts, free_errors, npts, supported_errors)
hold on
legend('clamped error','free error','supported plate')
xlabel('N')
ylabel('Relative error')
title('Analytic solution test')