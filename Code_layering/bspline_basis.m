function basis = bspline_basis(order,interval,knots)

m = order; %put your desired order of B-splines, e.g. m=3
sp = [0 1];%interval
L = interval;  %number of discrete point in interval
x = linspace(sp(1),sp(2),L); %discrete variable
b = 1/(knots - order + 1);  %distance between two adjacent knots
t = sp(1):b:sp(2); %knot sequence

tleft(1:m-1) = sp(1)-(m-1)*b:b:sp(1)-b;  %left knots for ESEP type basis
tright(1:m-1) = sp(2)+b:b:sp(2)+(m-1)*b; %right knots for ESEP type basis
t = [tleft t tright]; %extended knot sequence
basis = tspline(x,t,m);


