function [r, d, d2] = heart(t)

x = 16*sin(t).^3;
y = 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t);

x1 = 48*sin(t).^2.*cos(t);
y1 = -13*sin(t) + 10*sin(2*t) + 6*sin(3*t) + 4*sin(4*t);

x2 = 96*sin(t).*cos(t).^2 - 48*sin(t).^3;
y2 = -13*cos(t) + 20*cos(2*t) + 18*cos(3*t) + 16*cos(4*t);

r = [(x(:)).'; (y(:)).'];
d = [(x1(:)).'; (y1(:)).'];
d2 = [(x2(:)).'; (y2(:)).'];

end