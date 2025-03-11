 




% define the geometry of the circle 
function [r, d, d2]= circle(t,R) 

x = R*cos(t);
y = R*sin(t);

dx = -R*sin(t);
dy = R*cos(t);

dxx = R*cos(t);
dyy = R*sin(t);

r = [(x(:)).';(y(:)).'];
d = [(dx(:)).'; (dy(:)).'];
d2 = [(dxx(:)).'; (dyy(:)).'];


end