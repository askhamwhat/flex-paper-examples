
h = 0.01;

[tx1,tx2] = meshgrid(-4*h:h:4*h, 0:h:8*h);
x = tx1 + 1;
y = tx2 + 3;

val = x.^3+y.^3+x.*y;

d1f = [-49/20	6	-15/2	20/3	-15/4	6/5	-1/6 0 0] .';
d2f = [469/90	-223/10	879/20	-949/18	41	-201/10	1019/180	-7/10 0].';
d3f = [-801/80	349/6	-18353/120	2391/10	-1457/6	4891/30	-561/8	527/30	-469/240].';

d2c = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560];

d1dy1 = zeros(9);
d1dy1(:,5) = d1f / h;

d2dy2 = zeros(9);
d2dy2(:,5) = d2f / h^2;

d3dy3 = zeros(9);
d3dy3(:,5) = d3f / h^3;

d2dx2 = zeros(9);
d2dx2(1,:) = d2c / h^2;

err1 = 28 - sum(val.*d1dy1,'all')
err2 = 18 - sum(val.*d2dy2,'all')
err3 = 6 - sum(val.*d3dy3,'all')
err4 = 6 - sum(val.*d2dx2,'all')

val = x.^2.*y;

d3dx2dy = d1f*d2c / h^3;
err5 = 2 - sum(val.*d3dx2dy,'all')