function [r1,grad,hess,third] = planewave1(k,r,d)

r1 = exp(1i*sum(k.*r(:,:).*d(:)))';

gx =exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k);
gy =exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(2)*k);
grad = [gx; gy]';

hxx = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k).^2; 
hxy = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k).*(1i*d(2)*k); 
hyy = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(2)*k).^2;
hess = [hxx; hxy; hyy]';

txxx = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k).^3;
txxy = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k).^2*(1i*d(2)*k); 
txyy = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(1)*k)*(1i*d(2)*k).^2; 
tyyy = exp(1i*sum(k.*r(:,:).*d(:))).*(1i*d(2)*k).^3; 
third = [txxx; txxy; txyy; tyyy]';

end