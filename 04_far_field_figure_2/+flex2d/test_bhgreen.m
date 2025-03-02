%test_hkdiffgreen test the difference kernel H0 - K0 functions 
% using finite differences and high precision values
% 

clearvars; clc
rng(12345);

ns = 2;
nt = 4;

src0 = randn(2,ns)*1e-1;
targ0 = randn(2,nt)*1e-1;

% test derivatives

[val0,grad0,hess0,der30,der40,der50] = chnk.flex2d.bhgreen(src0,targ0);

for j = 1:5
    h = 10^(-j);
    dx = h*[1;0];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41,der51] = chnk.flex2d.bhgreen(src0,targ1);

    errdx = norm(ones(size(val0)) - 2*(val1-val0)/h./(grad0(:,:,1)+grad1(:,:,1)));
    errdxx = norm(ones(size(val0)) - 2*(grad1(:,:,1)-grad0(:,:,1))/h./(hess0(:,:,1)+hess1(:,:,1)));
    errdxy = norm(ones(size(val0)) - 2*(grad1(:,:,2)-grad0(:,:,2))/h./(hess0(:,:,2)+hess1(:,:,2)));
    errdxxx = norm(ones(size(val0)) - 2*(hess1(:,:,1)-hess0(:,:,1))/h./(der30(:,:,1)+der31(:,:,1)));
    errdxxy = norm(ones(size(val0)) - 2*(hess1(:,:,2)-hess0(:,:,2))/h./(der30(:,:,2)+der31(:,:,2)));
    errdxyy = norm(ones(size(val0)) - 2*(hess1(:,:,3)-hess0(:,:,3))/h./(der30(:,:,3)+der31(:,:,3)));
    errdxxxx = norm(ones(size(val0)) - 2*(der31(:,:,1)-der30(:,:,1))/h./(der40(:,:,1)+der41(:,:,1)));
    errdxxxy = norm(ones(size(val0)) - 2*(der31(:,:,2)-der30(:,:,2))/h./(der40(:,:,2)+der41(:,:,2)));
    errdxxyy = norm(ones(size(val0)) - 2*(der31(:,:,3)-der30(:,:,3))/h./(der40(:,:,3)+der41(:,:,3)));
    errdxyyy = norm(ones(size(val0)) - 2*(der31(:,:,4)-der30(:,:,4))/h./(der40(:,:,4)+der41(:,:,4)));
    errdxxxxx = norm(ones(size(val0)) - 2*(der41(:,:,1)-der40(:,:,1))/h./(der50(:,:,1)+der51(:,:,1)));
    errdxxxxy = norm(ones(size(val0)) - 2*(der41(:,:,2)-der40(:,:,2))/h./(der50(:,:,2)+der51(:,:,2)));
    errdxxxyy = norm(ones(size(val0)) - 2*(der41(:,:,3)-der40(:,:,3))/h./(der50(:,:,3)+der51(:,:,3)));
    errdxxyyy = norm(ones(size(val0)) - 2*(der41(:,:,4)-der40(:,:,4))/h./(der50(:,:,4)+der51(:,:,4)));
    errdxyyyy = norm(ones(size(val0)) - 2*(der41(:,:,5)-der40(:,:,5))/h./(der50(:,:,5)+der51(:,:,5)));

    dx = h*[0;1];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41,der51] = chnk.flex2d.bhgreen(src0,targ1);

    errdy = norm(ones(size(val0)) - 2*(val1-val0)/h./(grad0(:,:,2)+grad1(:,:,2)));
    errdyy = norm(ones(size(val0)) - 2*(grad1(:,:,2)-grad0(:,:,2))/h./(hess0(:,:,3)+hess1(:,:,3)));
    errdyyy = norm(ones(size(val0)) - 2*(hess1(:,:,3)-hess0(:,:,3))/h./(der30(:,:,4)+der31(:,:,4)));
    errdyyyy = norm(ones(size(val0)) - 2*(der31(:,:,4)-der30(:,:,4))/h./(der40(:,:,5)+der41(:,:,5)));
    errdyyyyy = norm(ones(size(val0)) - 2*(der41(:,:,5)-der40(:,:,5))/h./(der50(:,:,6)+der51(:,:,6)));

    fprintf('%5.2e : err in dx\n',errdx)    
    fprintf('%5.2e : err in dy\n',errdy)    
    fprintf('%5.2e : err in dxx\n',errdxx)    
    fprintf('%5.2e : err in dxy\n',errdxy)    
    fprintf('%5.2e : err in dyy\n',errdyy)    
    fprintf('%5.2e : err in dxxx\n',errdxxx)    
    fprintf('%5.2e : err in dxxy\n',errdxxy)    
    fprintf('%5.2e : err in dxyy\n',errdxyy)    
    fprintf('%5.2e : err in dyyy\n',errdyyy)    
    fprintf('%5.2e : err in dxxxx\n',errdxxxx)    
    fprintf('%5.2e : err in dxxxy\n',errdxxxy)    
    fprintf('%5.2e : err in dxxyy\n',errdxxyy)    
    fprintf('%5.2e : err in dxyyy\n',errdxyyy)    
    fprintf('%5.2e : err in dyyyy\n',errdyyyy)    
    fprintf('%5.2e : err in dxxxxx\n',errdxxxxx)    
    fprintf('%5.2e : err in dxxxxy\n',errdxxxxy)    
    fprintf('%5.2e : err in dxxxyy\n',errdxxxyy)    
    fprintf('%5.2e : err in dxxyyy\n',errdxxyyy)    
    fprintf('%5.2e : err in dxyyyy\n',errdxyyyy)    
    fprintf('%5.2e : err in dyyyyy\n',errdyyyyy)    
end

