function A = interset_area(Xg,Yg,lat,lon),

xsq = [Xg-0.5; Xg-0.5; Xg+0.5; Xg+0.5; Xg-0.5;];
ysq = [Yg-0.5; Yg+0.5; Yg+0.5; Yg-0.5; Yg-0.5;];

if ~ispolycw(lat, lon)
    [lat, lon] = poly2cw(lat, lon);
end;
if ~ispolycw(xsq, ysq)
    [xsq, ysq] = poly2cw(xsq, ysq);
end;

[xb, yb] = polybool('intersection', xsq, ysq, lat, lon);
A        = polyarea(xb(~isnan(xb)),yb(~isnan(yb)));
