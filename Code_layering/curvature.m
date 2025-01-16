function [xnormals,ynormals,curvaturemap,surprisalmap,angelmap] = curvature(x, y, signed, resolution)

x = x(:); y = y(:);

if length(x) ~= length(y)
error('the x and y vectors defining the contour must have equal lengths!');
end;

N = length(x);
b = 1; % spread term for the Von Mises distribution

if (x(1) == x(N) & y(1) == y(N))
N = N-1;
open = 0;
else
open = 1;
end

normalmap = [];
angelmap = [];
surprisalmap = [];
curvaturemap = [];
for j = 1:N

% define the previous and next points relative to the current one
xprev = 0;
yprev = 0;
xnext = 0;
ynext = 0;

for k = 1:resolution
prev(k) = j-k;
next(k) = j+k;
if(prev(k) < 1) prev(k) = prev(k) + N; end
if(next(k) > N) next(k) = next(k) - N; end

xprev = xprev + x(prev(k));
yprev = yprev + y(prev(k));
xnext = xnext + x(next(k));
ynext = ynext + y(next(k));
end

xprev = xprev/k;
yprev = yprev/k;
xnext = xnext/k;
ynext = ynext/k;


% vector _from_ the previous point; vector _to_ the next point
vecprev = [x(j), y(j)] - [xprev, yprev];
vecnext = [xnext, ynext] - [x(j), y(j)];

% compute the magnitude of the turning angle using dot product
alpha = acos( dot(vecprev, vecnext) / (norm(vecprev)*norm(vecnext)) );

% compute the sign of turning using cross product
% (assumes a counterclockwise sampling of the contour)
cp = cross( [vecprev, 0], [vecnext,0] );
alpha = sign(cp(3))*alpha;

% compute the surprisal using -log von mises
if(signed)
surprisal = -log( exp(b*cos(alpha - (2*pi/(N/resolution))) )/( 2*pi*besseli(0,b) ) );
else
surprisal = -log( exp(b*cos(alpha) )/( 2*pi*besseli(0,b) ) );
end;

% turning angles not defined near the end points of an open curve
if(open & (j - resolution < 1 | j + resolution > N)) surprisal = 0; end

% compute the tangent and normal vectors
tangvec = vecprev + vecnext;
tangvec = tangvec/norm(tangvec);
normvec = [[0 1; -1 0]*tangvec']';

normalmap = [normalmap; normvec];
surprisalmap = [surprisalmap; surprisal];
curvaturemap = [curvaturemap; alpha];
u  = [1 0 0];
v  = [normvec 0 ];
ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
angelmap = [angelmap, ThetaInDegrees];
end;

% scale the size of histogram bars
surprisalmapS = (1/range(surprisalmap(k+1:end-k)))*(surprisalmap(k+1:end-k) - min(surprisalmap(k+1:end-k)));
surprisalmapS = (surprisalmapS+0.2)./max(surprisalmapS+0.2);

surprisalmap = [ones(k,1)*surprisalmapS(1);surprisalmapS;ones(k,1)*surprisalmapS(end)];
curvaturemap = [ones(k,1)*curvaturemap(k+1);curvaturemap(k+1:end-k);ones(k,1)*curvaturemap(end-k)];
needlesize = max(range(x), range(y))/1;

% define the histogram normal to the shape
% normalmap(:,1) = surprisalmap.*normalmap(:,1);
% normalmap(:,2) = surprisalmap.*normalmap(:,2);
normalmap0 = normalmap; 
normalmap(:,1) = curvaturemap.*normalmap(:,1);
normalmap(:,2) = curvaturemap.*normalmap(:,2);
% xnormals = [x(1:N) - needlesize*normalmap(:,1), x(1:N) + needlesize*normalmap(:,1)];
% ynormals = [y(1:N) - needlesize*normalmap(:,2), y(1:N) + needlesize*normalmap(:,2)];
xnormals = [x(1:N) , x(1:N) + needlesize*normalmap(:,1)];
ynormals = [y(1:N) , y(1:N) + needlesize*normalmap(:,2)];
% plot
% figure, plot(x,y,'r'), axis equal;
% hold on, plot(xnormals', ynormals', 'k-');

% 
% % scale the size of histogram bars
% surprisalmap = (1/range(surprisalmap))*(surprisalmap(:) - min(surprisalmap));
% needlesize = max(range(x), range(y))/10;
% 
% % define the histogram normal to the shape
% normalmap(:,1) = surprisalmap.*normalmap(:,1);
% normalmap(:,2) = surprisalmap.*normalmap(:,2);
% xnormals = [x(1:N), x(1:N) + needlesize*normalmap(:,1)];
% ynormals = [y(1:N), y(1:N) + needlesize*normalmap(:,2)];
% 
% % plot
% figure, plot(x,y,'r'), axis equal;
% hold on, plot(xnormals', ynormals', 'k-');