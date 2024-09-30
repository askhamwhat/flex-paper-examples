% % Plane wave info

zk = 1.4;               % our k (wave number)
nu = 1/3;
theta = pi/4;
d = -[cos(theta) sin(theta)];

% Generating chunkr objects

valid_positions = false;
box_l = 0; box_r = 49;
box_d = 0; box_u = 49;
s = 10;
[xs, ys] = meshgrid(box_l+s/2:s:box_r-s/2,box_d+s/2:s:box_u-s/2);
centers = [xs(:) ys(:)];
centers = centers + (rand(size(centers))-0.5)*s/2;
centers(1,:)  = [18,18];
n_objs = numel(xs(:));

% chnkrs = createArray(n_objs,1,'chunker');
chnkrs = [];
cparams.maxchunklen = 100;

% Plotting and chunkr creation done in the same for loop


clf
figure(1) 
for i = 1:n_objs
    narms = randi([3,6]);
    narms = 4;
    chnkr = chunkerfunc(@(t) starfish(t,narms), cparams);
    scale = rand/2+1;
    scale = 1.0885;
    rot = rand*2*pi;
    rot = 1.5457;
    chnkr = chnkr.move([0;0],centers(i,:)',rot,scale);
    chnkrs = [chnkrs,chnkr];

    plot(chnkr, '-x'); hold on
    quiver(chnkr); hold on
end
axis equal
drawnow

% Defining kernels

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
fkern1_trans = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate hilbert unsubtract', coefs);  % K11 without the hilbert subtraction
fkern3_trans = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate K21 first part unsubtract', coefs);             % singularity subtration kernel in K21 (including swapping its Asmyptotics expansions)

%Building lefthand side 

start = tic;

lhs = [];

for j = 1:n_objs

    src_chnkr = chnkrs(j);
    opts2 = [];
    opts2.sing = 'pv';
    H = chunkermat(src_chnkr, hilbert, opts2); 
    wts = src_chnkr.wts(:)';

    col = [];

    for i = 1:n_objs

        targ_chnkr = chnkrs(i);

        if i == j % build scattering matrix

            sysmat = chunkermat(targ_chnkr,fkern, opts);
            sysmat1 = chunkermat(targ_chnkr, fkern1, opts);
            sysmat2 = chunkermat(targ_chnkr, fkern2, opts);

            K21 = chunkermat(targ_chnkr, fkern3, opts);
            K21second = chunkermat(targ_chnkr, fkern4, opts);
            K21hilbert = chunkermat(targ_chnkr, fkern5, opts);
            K22second = chunkermat(targ_chnkr, fkern6, opts);
            D = chunkermat(targ_chnkr, double, opts);                                            % Assemble hilbert transforms

            kappa = signed_curvature(targ_chnkr);
            kappa = kappa(:)';

            hilb = sysmat1*H - ((1+nu)/2).*(D*D)- ((1+nu)*nu/2).*(D*D);
            hilb2 = sysmat2*H ;

            mat1 =  sysmat(1:2:end, 1:2:end);
            mat4 = sysmat(2:2:end, 2:2:end);

            sysmat(1:2:end, 1:2:end) = mat1 + hilb;
            sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa.').*(K21hilbert*H + K21second);
            sysmat(2:2:end, 2:2:end) = mat4 + (kappa.').*(K22second);

            A = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];                                     % jump matrix (for exterior problem)
            M = kron(eye(targ_chnkr.npt), A);
            Sinv =  M + sysmat ;

            col = [col; Sinv];

            disp(['S_',num2str(i)])

        else % build translation matrix 

            sysmat = fkern(src_chnkr,targ_chnkr);
            sysmat(:,1:2:end) = sysmat(:,1:2:end).*wts;
            sysmat(:,2:2:end) = sysmat(:,2:2:end).*wts;
            sysmat1 = fkern1_trans(src_chnkr,targ_chnkr).*wts;
            sysmat2 = fkern2(src_chnkr,targ_chnkr).*wts;

            K21 = fkern3_trans(src_chnkr,targ_chnkr).*wts;
            K21second = fkern4(src_chnkr,targ_chnkr).*wts;
            K21hilbert = fkern5(src_chnkr,targ_chnkr).*wts;
            K22second = fkern6(src_chnkr,targ_chnkr).*wts;


            kappa = signed_curvature(targ_chnkr);
            kappa = kappa(:)';

            hilb = sysmat1*H ;
            hilb2 = sysmat2*H ;

            mat1 =  sysmat(1:2:end, 1:2:end);
            mat4 = sysmat(2:2:end, 2:2:end);

            sysmat(1:2:end, 1:2:end) = mat1 + hilb;
            sysmat(2:2:end, 1:2:end) = K21 +  hilb2 + (kappa.').*(K21hilbert*H + K21second);
            sysmat(2:2:end, 2:2:end) = mat4 + (kappa.').*(K22second);

            T =  sysmat;

            col = [col; T];

            disp(['T_',num2str(i),',',num2str(j)])

        end

    end

    lhs = [lhs col];
end 

t1 = toc(start);
fprintf('%5.2f s : time to assemble matrix\n',t1)

% Evaluating righthand side(s)

start = tic;
rhs = [];

src = [18.1;17.9];
for i = 1:n_objs

    chnkr = chnkrs(i);

    [~, ~, hess, third] = planewave1(zk, chnkr.r(:,:), d);
    [~,~,hess,third,~] = chnk.flex2d.helmdiffgreen(zk,src,chnkr.r(:,:));
    
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
    
    data = zeros(2*chnkr.npt, 1); 
    data(1:2:2*chnkr.npt) = -firstbc ; 
    data(2:2:2*chnkr.npt) = -secondbc;
    
    rhs = [rhs; data];

end

t1 = toc(start);
fprintf('%5.2f s : time to evaluate BCs \n',t1) 


start = tic;

sol = lhs\rhs;

t1 = toc(start);
fprintf('%5.2f s : time to solve linear system \n',t1)  

start1 = tic;

% Plotting the solution

xs = box_l-s/2:0.25:box_r+s/2;                                    % generate some targets
ys = box_d-s/2:0.25:box_u+s/2;
xs =linspace(14,22,100);
ys = xs;
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);

in = 0;
for i = 1:n_objs
    in = in + chunkerinterior(chnkrs(i), targets); 
end
out = ~in; 

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

uscat = 0;
count = 0;
for i = 1:n_objs
    rho1 = sol(count+1:2:count+2*chnkrs(i).npt); 
    rho2 = sol(count+2:2:count+2*chnkrs(i).npt); 
    count = count + 2*chnkrs(i).npt;

    H = chunkermat(chnkrs(i), hilbert, opts2);                                              % Assemble hilbert transforms
    
    uscat =uscat + chunkerkerneval(chnkrs(i), ikern1,rho1, targets(:, out)) +...
        chunkerkerneval(chnkrs(i), ikern3,H*rho1, targets(:, out)) + ...
        chunkerkerneval(chnkrs(i), ikern2, rho2, targets(:,out));

    disp(['Object ', num2str(i), 'kernel evaluation complete'])
end 
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

 [uinc,~,~,~] = chnk.flex2d.helmdiffgreen(zk,src,targets(:,out));
%[uinc,~, ~] = planewave1(zk,targets(:,out),d);


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
set(h,'EdgeColor','None'); hold on;
title("Incident field", 'FontSize',16)
plot(chnkrs,'w-','LineWidth',2);
clim([-1,1]);
colorbar
set(gca,'fontsize',16)
% set(gca,'TickLabelInterpreter','latex')

nexttile
h = pcolor(X,Y,real(us));
set(h,'EdgeColor','None'); hold on;
title("Scattered field", 'FontSize',16)
plot(chnkrs,'w-','LineWidth',2);
clim([-1,1]);
colorbar
set(gca,'fontsize',16)
% set(gca,'TickLabelInterpreter','latex')

nexttile
h = pcolor(X,Y,log10(abs(utot)));
set(h,'EdgeColor','None'); hold on;
title("Total field", 'FontSize',16)
plot(chnkrs,'w-','LineWidth',2);
clim([-1.3,1.3]);
colorbar
set(gca,'fontsize',16)
% set(gca,'TickLabelInterpreter','latex')

% us = NaN*zeros(na, 1);
% ui = NaN*zeros(na,1);
% us(out) = uscat;
% ui(out) = uinc;
% us = reshape(us,[numel(xs) numel(ys)]);
% ui = reshape(ui,[numel(xs) numel(ys)]);
% utot = us + ui;
% 
% 
% figure(4)
% tiledlayout(1,3)
% nexttile
% h = pcolor(X,Y,real(ui));
% set(h,'EdgeColor','None'); hold on;
% title("Incident field", 'FontSize',16)
% plot(chnkrs,'k-','LineWidth',1);
% clim([-1,1]);
% colorbar
% set(gca,'fontsize',16)
% % set(gca,'TickLabelInterpreter','latex')
% 
% nexttile
% h = pcolor(X,Y,real(us));
% set(h,'EdgeColor','None'); hold on;
% title("Scattered field", 'FontSize',16)
% plot(chnkrs,'k-','LineWidth',1);
% clim([-1,1]);
% colorbar
% set(gca,'fontsize',16)
% % set(gca,'TickLabelInterpreter','latex')
% 
% nexttile
% h = pcolor(X,Y,real(utot));
% set(h,'EdgeColor','None'); hold on;
% title("Total field", 'FontSize',16)
% plot(chnkrs,'k-','LineWidth',1);
% clim([-1.3,1.3]);
% colorbar
% set(gca,'fontsize',16)

