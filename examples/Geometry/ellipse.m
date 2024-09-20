function [r, d, d2] = ellipse(t)

c = 2;

x = c.*cos(t);
y = sin(t);

x1 = -c.*sin(t);
y1 = cos(t);

x2 = -c.*cos(t);
y2 = -sin(t);

r = [(x(:)).'; (y(:)).'];
d = [(x1(:)).'; (y1(:)).'];
d2 = [(x2(:)).'; (y2(:)).'];

end