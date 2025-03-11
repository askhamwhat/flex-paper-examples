function [r, d, d2] = blob(t)

a = 2;
b = 1;
c = -0.4;

x = a*cos(t);
y = b*sin(t) + c*cos(t).^2;

x1 = -a*sin(t);
y1 = b*cos(t) - 2*c*cos(t).*sin(t);

x2 = -a*cos(t);
y2 = -b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;

r = [(x(:)).'; (y(:)).'];
d = [(x1(:)).'; (y1(:)).'];
d2 = [(x2(:)).'; (y2(:)).'];

end