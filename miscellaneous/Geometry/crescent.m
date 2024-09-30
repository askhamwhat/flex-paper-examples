function [r, d, d2] = crescent(t)

a = -0;
b = 0.5;
c = -0.2;

x = cos(t) + a*sin(t).^2;
y = b*sin(t) + c*cos(t).^2;

x1 = -sin(t) + 2*a*sin(t).*cos(t);
y1 = b*cos(t) - 2*c*cos(t).*sin(t);

x2 = - cos(t) + 2*a*cos(t).^2 - 2*a*sin(t).^2;
y2 = - b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;

r = 2*[(x(:)).'; (y(:)).'];
d = 2*[(x1(:)).'; (y1(:)).'];
d2 = 2*[(x2(:)).'; (y2(:)).'];

end