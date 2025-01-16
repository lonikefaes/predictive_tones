function varargout = clipPointsSphare3d(points,P,radius)
%CLIPPOINTS3D Clip a set of points by a box
%
%   CLIP = clipPoints3d(POINTS, BOX);
%   Returns the set of points which are located inside of the BOX.
%
%   [CLIP IND] = clipPoints2d(POINTS, BOX);
%   Also returns the indices of clipped points.
%
%   See also
%   points3d, boxes3d
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.

% get bounding box limits


% compute indices of points inside visible area
ind = ((points(:,1)- P(1,1)).^2 + (points(:,2)-P(1,2)).^2 + (points(:,3)-P(1,3)).^2) <= radius.^2;

% keep only points inside box
ind    = find(ind==1);
points = points(ind, :);

% process output arguments
varargout{1} = points;
if nargout == 2
    varargout{2} = ind;
end
