function tf = ispolycw(x, y)
%ISPOLYCW True if polygon vertices are in clockwise order
%
%   TF = ISPOLYCW(X, Y) returns true if the polygonal contour vertices 
%   represented by X and Y are ordered in the clockwise direction.  X and Y
%   are numeric vectors with the same number of elements.
%
%   Alternatively, X and Y can contain multiple contours, either in
%   NaN-separated vector form or in cell array form.  In that case,
%   ISPOLYCW returns a logical array containing one true or false value
%   per contour.
%
%   ISPOLYCW always returns true for polygonal contours containing two or
%   fewer vertices.
%
%   Vertex ordering is not well defined for self-intersecting polygonal
%   contours.  For such contours, ISPOLYCW returns a result based on the
%   order or vertices immediately before and after the left-most of the 
%   lowest vertices.  In other words, of the vertices with the lowest Y
%   value, find the vertex with the lowest X value.  For a few special
%   cases of self-intersecting contours, the vertex ordering cannot be
%   determined using only the left-most of the lowest vertices; for these
%   cases, ISPOLYCW uses a signed area test to determine the ordering.
%
%   Class Support
%   -------------
%   X and Y may be any numeric class.
%
%   Example
%   -------
%   Orientation of a square:
%
%       x = [0 1 1 0 0];
%       y = [0 0 1 1 0];
%       ispolycw(x, y)                     % Returns 0
%       ispolycw(fliplr(x), fliplr(y))     % Returns 1
%
%   See also POLY2CW, POLY2CCW, POLYBOOL.

% Copyright 2004-2009 The MathWorks, Inc.
% $Revision: 1.1.4.3 $  $Date: 2009/08/11 15:44:28 $

if isempty(x)
   tf = true;
   return;
end

if ~iscell(x)
   checkxy(x, y, mfilename, 'X', 'Y', 1, 2)
   is_row = (size(x,1) == 1);
   [x, y] = polysplit(x, y);
   if is_row
      x = x';
      y = y';
   end
end

tf = false(size(x));
for k = 1:numel(x)
   tf(k) = isContourClockwise(x{k}, y{k});
end

%----------------------------------------------------------------------
function tf = isContourClockwise(x, y)

if numel(x) <= 1
    tf = true;
    return;
end

is_closed = (x(1) == x(end)) && (y(1) == y(end));
if is_closed
    x(end) = [];
    y(end) = [];
end

[x, y] = removeDuplicates(x, y);
num_vertices = numel(x);
if num_vertices <= 2
    tf = true;
    return;
end

idx = findExtremeVertices(x, y);

if numel(idx) > 1
    % The same extreme vertex appears multiple, nonsuccessive times in
    % the vertex list.  Use signed area test.
    tf = signedArea(x, y) <= 0;
    return;
end

% Find the three vertices we are interested in: the left-most of the
% lowest vertices, as well as the ones immediately before and after it.
p = mod((idx - 1) + [-1, 0, 1], num_vertices) + 1;
xx = x(p);
yy = y(p);

if ~isfloat(xx)
    xx = double(xx);
end
if ~isfloat(yy)
    yy = double(yy);
end

ux = xx(2) - xx(1);
uy = yy(2) - yy(1);

vx = xx(3) - xx(2);
vy = yy(3) - yy(2);

a = ux*vy;
b = uy*vx;
if a == b
    % The left-most lowest vertex is the end-point of a kind of linear
    % "spur."  The contour doubles back on itself, such as in this case:
    % x = [0 1 1 0 0 -1 0];
    % y = [0 0 1 1 0 -1 0];
    % The left-most lowest vertex is (-1,-1), but we since this vertex
    % is the end-point of a spur, we can't tell the direction from it.
    % Use the signed polygon test.
    tf = signedArea(x, y) <= 0;
else
    tf = a < b;
end

%----------------------------------------------------------------------
function [xout, yout] = removeDuplicates(x, y)
num_vertices = numel(x);
k1 = [2:num_vertices 1];
k2 = 1:num_vertices;
dups = (x(k1) == x(k2)) & (y(k1) == y(k2));
xout = x;
yout = y;
xout(dups) = [];
yout(dups) = [];

%----------------------------------------------------------------------
function idx = findExtremeVertices(x, y)
% Return the indices of all the left-most lowest vertices in (x,y).

% Find the vertices with the minimum y.
idx = find(y == min(y));

x_subset = x(idx);
idx2 = (x_subset == min(x_subset));

idx = idx(idx2);

%----------------------------------------------------------------------
function a = signedArea(x, y)
% a = signedArea(x,y) returns twice the signed area of the polygonal
% contour represented by vectors x and y.  Assumes (x,y) is NOT closed.

% Reference: 
% http://geometryalgorithms.com/Archive/algorithm_0101/algorithm_0101.htm

x = x - mean(x);
n = numel(x);
if n <= 2
    a = 0;
else
    i = [2:n 1];
    j = [3:n 1 2];
    k = (1:n);
    a = sum(x(i) .* (y(j) - y(k)));
end
