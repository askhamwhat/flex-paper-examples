%test_hkdiffgreen test the difference kernel H0 - K0 functions 
% using finite differences and high precision values
% 

clearvars; clc
rng(12345);

ns = 2;
nt = 4;
k = 1.3;

src0 = randn(2,ns)*1e-1;
targ0 = randn(2,nt)*1e-1;

% test derivatives

ifr2logr = false;
[val0,grad0,hess0,der30,der40,der50] = flex2d.hkdiffgreen(k,src0,targ0,ifr2logr);

for j = 1:5
    h = 10^(-j);
    dx = h*[1;0];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41,der51] = flex2d.hkdiffgreen(k,src0,targ1,ifr2logr);

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
    [val1,grad1,hess1,der31,der41,der51] = flex2d.hkdiffgreen(k,src0,targ1,ifr2logr);

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

%% some high precision calcs for comparison 

%% H_0(kr) - H_0(i kr)

% non power series eval 
ifr2logr = false;
srct = [0;0]; targt = [2.3;1.7];
kt = 1.1; 
valt = -0.086361033112625579806-1i*0.076378624373238939605;
gradxt = 0.084324892535124417680-0.062550394133266982860*1i;
hessxyt = 0.009494827062481582523+0.070174358169356783483*1i;
der3xxyt = -0.049828271176969417642-0.012236588644907110535*1i;
der4xxyyt = 0.0097294952080452932969-0.0082189736022730294709*1i;
der5xxxyyt = -0.0093420109741399530893-0.0099716715832952414166*1i;

[val,grad,hess,der3,der4,der5] = flex2d.hkdiffgreen(kt,srct,targt,ifr2logr);

abs(val-valt)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3xxyt)/abs(der3xxyt)
abs(der4(3)-der4xxyyt)/abs(der4xxyyt)
abs(der5(3)-der5xxxyyt)/abs(der5xxxyyt)

%% power series eval 

ifr2logr = false;
srct = [0;0]; targt = [2.3;1.7]*10^(-5);
kt = 1.1; 
valt = -9.0442561921*10^-10 + 0.24999999993813875000*1i;
gradxt = -0.000048645479362495222741-3.478749999569600353*10^-6*1i;
hessxyt = 0.092051094286316139716+1.7889471874*(10^-11)*1i;
der3xxyt = -1174.24591414456647061+7.7780312485*(10^-7)*1i;
der4xxyyt = 1.9489279020645010424*10^8+0.04575312499*1i;
der5xxxyyt = -2.1938397589649620042*10^13-0.*10^-7*1i;

[val,grad,hess,der3,der4,der5] = flex2d.hkdiffgreen(kt,srct,targt,ifr2logr);

abs(val-valt)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt) 
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3xxyt)/abs(der3xxyt)
abs(der4(3)-der4xxyyt)/abs(der4xxyyt) 
abs(der5(3)-der5xxxyyt)/abs(der5xxxyyt)

%% H_0(kr) - H_0(i kr) - r2 log (r)

%% power series eval 

ifr2logr = true;
srct = [0;0]; targt = 10^(-5)*[2.3;1.7];
kt = 1.1;
[val,grad,hess,der3,der4,der5] = flex2d.hkdiffgreen(kt,srct,targt,ifr2logr);

valt = -8.038841260*10^-11+0.24999999993813875000*1i;
gradxt = -4.5206197793280228361*10^-6-3.4787499995696003531*10^-6*1i;
hessxyt = -2.23557343249081186853*10^-20+1.78894718735244465440091814573*10^-11*1i;
der3xxyt = -2.17633187153*10^-15+7.7780312485286816821*10^-7*1i;
der4xxyyt = -1.64526768959*10^-10+0.045753124988678617985*1i;
der5xxxyyt = -8.9028273044542249271*10^-6-6.366547342908069210*10^-7*1i;

abs(val-valt)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3xxyt)/abs(der3xxyt)
abs(der4(3)-der4xxyyt)/abs(der4xxyyt)
abs(der5(3)-der5xxxyyt)/abs(der5xxxyyt)

%% non power series 

ifr2logr = true;
srct = [0;0]; targt = [2.3;1.7];
kt = 1;
[val,grad,hess,der3,der4,der5] = flex2d.hkdiffgreen(kt,srct,targt,ifr2logr);

der5xxxyyt = 0.0065718665239521022295-0.0088604567871282522532*1i;

abs(der5(3)-der5xxxyyt)/abs(der5xxxyyt)
