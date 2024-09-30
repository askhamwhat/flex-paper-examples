function [r, d, d2] = rotateline(t, theta)

x = cos(theta).*t;
y = sin(theta).*t;

dx = cos(theta).*ones(size(x));
dy = sin(theta).*ones(size(x));

dxx = 0.*t;
dyy = 0.*t;

r = [(x(:)).' ; (y(:)).'];
d = [(dx(:)).' ; (dy(:)).'];
d2 = [(dxx(:)).' ; (dyy(:)).'];

end