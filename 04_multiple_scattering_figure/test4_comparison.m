% Parameters

n_objs = 2; % number of objects
source_loc = [-3;0]; % location of test source

zk = 1;  % our k (wave number)
nu = 1/3; % Poisson ratio
cparams = [];
cparams.maxchunklen = 0.5 / zk; % max chunk size
narms = 5; % number of arms on the starfish

s = 5; % spacing
box_size = s*ceil(sqrt(n_objs))-s;
[xs, ys] = meshgrid(0:s:box_size,0:s:box_size);
centers = [xs(:) ys(:)];

temp = [centers, abs(centers(:,1) - centers(:,2))];
temp = sortrows(temp, 3);
centers = temp(1:n_objs,1:2);

% Generating chnkr objects

%chnkrs = createArray(n_objs,1,'chunker');
chnkrs = [];
for i = 1:n_objs
    chnkr = chunkerfunc(@(t) starfish(t,narms), cparams);
    chnkr = chnkr.move([0;0],centers(i,:)',2*i,1);
    %chnkrs(i) = chnkr;
    chnkrs = [chnkrs,chnkr];
end

clf
figure(1)
plot(chnkrs, '-x'); hold on
quiver(chnkrs); hold on
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

% Building lefthand side 

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

for i = 1:n_objs

    chnkr = chnkrs(i);

    [~, ~, hess, third, ~] = chnk.flex2d.helmdiffgreen(zk, source_loc, chnkr.r);
    [~, ~, hessK, thirdK, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), source_loc, chnkr.r);
    
    nx = chnkr.n(1,:).'; 
    ny = chnkr.n(2,:).';
    
    dx = chnkr.d(1,:).';
    dy = chnkr.d(2,:).';
    
    ds = sqrt(dx.*dx+dy.*dy);
    taux = (dx./ds);                                                                       % normalization
    tauy = (dy./ds);
    
    kappa = signed_curvature(chnkr);
    kappa = kappa(:)';
    
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
            (1-coefs(1))*diag(kappa)*(1/(2*zk^2).*(hess(:, :, 1).*taux.*taux + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*tauy.*tauy)-...
             1/(2*zk^2).*(hessK(:, :, 1).*taux.*taux + hessK(:, :, 2).*(2*taux.*tauy) + ...
            hessK(:, :, 3).*tauy.*tauy)-...
            (1/(2*zk^2).*(hess(:, :, 1).*nx.*nx + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*ny.*ny)-...
           1/(2*zk^2).*(hessK(:, :, 1).*nx.*nx + hessK(:, :, 2).*(2*nx.*ny) + hessK(:, :, 3).*ny.*ny)));

    data = zeros(2*chnkr.npt, 1); 
    data(1:2:2*chnkr.npt) = firstbc ; 
    data(2:2:2*chnkr.npt) = secondbc;
    
    rhs = [rhs; data];

end

t1 = toc(start);
fprintf('%5.2f s : time to evaluate BCs \n',t1) 

% Solving linear system

start = tic;

sol = lhs\rhs;

t1 = toc(start);
fprintf('%5.2f s : time to solve linear system \n',t1)  

start1 = tic;



% Plotting the solution

xs = -s/2:s/20:box_size+s/2;                                    % generate some targets
ys = -s/2:s/20:box_size+s/2; 
[X,Y] = meshgrid(xs, ys);
targets = [X(:).'; Y(:).'];
[~,na] = size(targets);
targets = [0;3];

start = tic;
in = 0;
for i = 1:n_objs
    in = in + chunkerinterior(chnkrs(i), targets); 
end
out = ~in; 
t1 = toc(start);
fprintf('%5.2f s : time to find interior points \n',t1)  

ikern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first', coefs);                              % build the kernel of evaluation          
ikern2 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
ikern3 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert',coefs);

%H = chunkermat(chnkr1, hilbert, opts2);                                              % Assemble hilbert transforms

uscat = 0;
count = 0;
start1 = tic;
for i = 1:n_objs
    rho1 = sol(count+1:2:count+2*chnkrs(i).npt); 
    rho2 = sol(count+2:2:count+2*chnkrs(i).npt); 
    count = count + 2*chnkrs(i).npt;

    H = chunkermat(chnkrs(i), hilbert, opts2);                                              % Assemble hilbert transforms
    
    uscat =uscat + chunkerkerneval(chnkrs(i), ikern1,rho1, targets(:, out)) +...
        chunkerkerneval(chnkrs(i), ikern3,H*rho1, targets(:, out)) + ...
        chunkerkerneval(chnkrs(i), ikern2, rho2, targets(:,out));

    disp(['Object ', num2str(i), ' kernel evaluation complete'])
end 


[val, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk, source_loc, targets(:,out));
[valK, ~, ~, ~, ~] = chnk.flex2d.helmdiffgreen(zk*(1i), source_loc, targets(:,out));

trueval = 1/(2*zk^2).*val - 1/(2*zk^2).*valK;

true_sol = zeros(na, 1);
utarg = zeros(na, 1);

utarg(out) = uscat;
true_sol(out) = trueval;

utarg = reshape(utarg,size(X));
true_sol = reshape(true_sol,size(X));

t2 = toc(start1);
fprintf('%5.2f s : time for kernel eval (for plotting)\n',t2)

figure(2)
tiledlayout(1,3)
nexttile
h = pcolor(X,Y,real(true_sol));
set(h,'EdgeColor','None'); hold on;
title("True solution")
plot(chnkrs,'w-','LineWidth',2);
colorbar

nexttile
h = pcolor(X,Y,real(utarg));
set(h,'EdgeColor','None'); hold on;
title("BIE solution")
plot(chnkrs,'w-','LineWidth',2);
colorbar

nexttile
uerr = utarg - true_sol;
uerr = reshape(uerr,size(X));
h = pcolor(X,Y,log10(abs(uerr) / max(abs(true_sol),[],'all')));
set(h,'EdgeColor','None'); hold on;
title("Log_{10} relative error (coupled)")
plot(chnkrs,'w-','LineWidth',2);
colorbar
