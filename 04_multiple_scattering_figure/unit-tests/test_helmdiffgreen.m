%test_helmdiffgreen test the difference kernel functions 
% using finite differences 
% 

clearvars; clc
rng(12345);

ns = 3;
nt = 4;
k = 1.3;

src0 = randn(2,ns);
targ0 = randn(2,nt);

% test derivatives

ifr2logr = false;
[val0,grad0,hess0,der30,der40,der50] = flex2d.helmdiffgreen(k,src0,targ0,ifr2logr);

for j = 1:5
    h = 10^(-j);
    dx = h*[1;0];
    targ1 = targ0 + dx;
    [val1,grad1,hess1,der31,der41,der51] = flex2d.helmdiffgreen(k,src0,targ1,ifr2logr);

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
    [val1,grad1,hess1,der31,der41,der51] = flex2d.helmdiffgreen(k,src0,targ1,ifr2logr);

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

%%
% some high precision calcs for comparison

clearvars
% from mathematica
srct = [0;0]; targt = 10^(-4)*[1;1];
kt = sqrt(2);
valt = -0.03670784159258519723+0.24999999750000000625 *1i;
gradxt = -0.00014535819351409084213-0.00002499999987500000021 *1i; 
hessxyt = 0.079577479012801035120+1.249999995833*1e-9*1i;
der3xxyt = 0.000070689659841738687819+0.000012499999916666666823*1i; 
der3xxxt = 1591.54965094568054413183407228+0.00003749999983333333359375*1i;
der4xxyyt = 7.95774786149136280529936105345*1e6+0.12499999875000000364583*1i;
[val,grad,hess,der3,der4] = flex2d.helmdiffgreen(kt,srct,targt);

abs(val-valt)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(1)-der3xxxt)/abs(der3xxxt)
abs(der4(3)-der4xxyyt)/abs(der4xxyyt)

%abs(der3(2)-der3xxyt)/abs(der3xxyt) % this derivative appears to be just hard to get
% condition number of previous times machine precision
%eps(1)*abs(targt(2)*der4xxyyt)/abs(der3xxyt)



%%
% high precision calcs for second difference 

clearvars
ifr2logr = true;
srct = [0;0]; targt = 10^(-4)*[1;1];
kt = sqrt(2);
valt = -0.036707827485462220021148457259+0.249999997500000006249999993056*1i;
der3xxxt = 0.000220026727186442996438555226113+0.000037499999833333333593749999792*1i;
der3xxyt = 0.000070689659841738687819+0.000012499999916666666823*1i;

[val,grad,hess,der3,der4] = flex2d.helmdiffgreen(kt,srct,targt,ifr2logr);
    
abs(val-valt)/abs(valt)
abs(der3(1)-der3xxxt)/abs(der3xxxt)
abs(der3(2)-der3xxyt)/abs(der3xxyt)

eps(1)*(abs(targt(2)*der4(2))+abs(targt(1)*der4(1)))/abs(der3xxxt)

srct = [0;0]; targt = 10^(-2)*[3;1];
kt = sqrt(2);
[val,grad,hess,der3,der4,der5] = flex2d.helmdiffgreen(kt,srct,targt,ifr2logr);
der3xxyt = 0.0024439996259246495294+0.0012494167265597222964*1i;

abs(der3(2)-der3xxyt)/abs(der3xxyt)
    
%
ifr2logr = true;
srct = [0;0]; targt = 10^(-6)*[3;1];
kt = sqrt(2);
[val,grad,hess,der3,der4] = flex2d.helmdiffgreen(kt,srct,targt,ifr2logr);
der3xxyt = 9.7749591459180136637*1e-7+1.2499999999941666667*1e-7*1i;
der4xxxxt = 2.8465440744964340700+0.3749999999971250000*1i;
der4xxxyt = -0.028647889759656676367-3.75000000*1e-13*1i;
hessxyt = 3.1473469169590242094*1e-12+3.749999999993750000*1e-13*1i;
valt = -0.03670782626160332761+0.24999999999875000000*1i;
gradxt = -3.6734135047637916559*1e-7 - 7.4999999999812500000*1e-7*1i;
abs(valt-val)/abs(valt)
abs(grad(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3xxyt)/abs(der3xxyt)
abs(der4(1)-der4xxxxt)/abs(der4xxxxt)
abs(der4(2)-der4xxxyt)/abs(der4xxxyt)

%% diff 


ifr2logr = false;
srct = [0;0]; targt = 10^(-6)*[3;1];
kt = sqrt(2);
[val,grad,hess,der3,der4,der5] = flex2d.helmdiffgreen(kt,srct,targt,ifr2logr);
[valK,gradK,hessK,der3K,der4K,der5K] = flex2d.helmdiffgreen(1i*kt,srct,targt,ifr2logr);

valt = -2.138020014*10^-11+1i*0.24999999999875000000;
gradxt = -0.0000123506552532202971795-1i*7.499999999981250000*(10^-7);
hessxyt = 0.095492965855137201461+1i*3.75000000*10^-13;
der3xxyt = -25464.79089470325372303145030186339119925+1i*1.2499999999941666666666726562*10^-7;
der4xxxyt = 1.1459155902616464175353399930964341754524*10^10-1i*3.74999999999250000*10^-13;
der5xxxyyt = -6.875493541569878505221990510956141724406*10^15-1i*3.74999999999156250*10^-7;

abs(val-valK-valt)/abs(valt)
abs(grad(1)-gradK(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessK(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3K(2)-der3xxyt)/abs(der3xxyt)
abs(der4(2)-der4K(2)-der4xxxyt)/abs(der4xxxyt)
abs(der5(3)-der5K(3)-der5xxxyyt)/abs(der5xxxyyt)

%% diff with diff

ifr2logr = true;
srct = [0;0]; targt = 10^(-6)*[3;1];
kt = sqrt(2);
[val,grad,hess,der3,der4,der5] = flex2d.helmdiffgreen(kt,srct,targt,ifr2logr);
[valK,gradK,hessK,der3K,der4K,der5K] = flex2d.helmdiffgreen(1i*kt,srct,targt,ifr2logr);

valt = -1.2244711683090424373835721417*10^-12+1i*0.2499999999987500000000015624999999991319; 
gradxt = -0.0000123506552532202971795-1i*7.499999999981250000*(10^-7);
hessxyt = 0.095492965855137201461+1i*3.75000000*10^-13;
der3xxyt = -1.004816226109327373898816282298*10^-17+1i*1.2499999999941666666666726562499999972222*10^-7;
der4xxxyt = -6.231031856692311957281356178160542931985*10^-12-1i*3.74999999999250000000000515624999999818*10^-13;
der5xxxyyt = -6.875493541569878505221990510956141724406*10^15-1i*3.74999999999156250*10^-7;

abs(val-valK-valt)/abs(valt)
abs(grad(1)-gradK(1)-gradxt)/abs(gradxt)
abs(hess(2)-hessK(2)-hessxyt)/abs(hessxyt)
abs(der3(2)-der3K(2)-der3xxyt)/abs(der3xxyt)
abs(der4(2)-der4K(2)-der4xxxyt)/abs(der4xxxyt)
abs(der5(3)-der5K(3)-der5xxxyyt)/abs(der5xxxyyt)