%MOVMAX   Moving maximum value.
%   Y = MOVMAX(X,K) for a vector X and positive integer scalar K computes a
%   centered moving maximum by sliding a window of length K along X. Each
%   element of Y is the local maximum of the corresponding values of X
%   inside the window, with Y the same size as X. When K is even, the
%   window is centered about the current and previous elements of X. The
%   sliding window is truncated at the endpoints where there are fewer than
%   K elements from X to fill the window.
%   
%   For N-D arrays, MOVMAX operates along the first array dimension whose
%   size does not equal 1.
%
%   Y = MOVMAX(X,[NB NF]) for a vector X and nonnegative integers NB and NF
%   computes a moving maximum along the length of X, returning the local
%   maximum of the previous NB elements, the current element, and the next
%   NF elements of X.
%
%   Y = MOVMAX(...,DIM) operates along dimension DIM of X.
%
%   MOVMAX(...,MISSING) specifies how NaN (Not-a-Number) values are treated
%   and can be one of the following:
%
%       'omitnan'      - (default) the maximum of any window containing NaN
%                        values is the maximum of all its non-NaN elements.
%                        If all elements are NaN, the result is NaN.
%       'includenan'   - the maximum of any window containing NaN values is
%                        also NaN.
%
%   MOVMAX(...,'Endpoints',ENDPT) controls how the maximum is calculated at
%   the endpoints of X, where there are not enough elements to fill the
%   window. ENDPT can be either a scalar numeric or logical value or one of
%   the following:
%
%       'shrink'    - (default) compute the maximum over the number of
%                     elements of X that are inside the window, effectively
%                     reducing the window size to fit X at the endpoints.
%       'fill'      - compute the maximum over the full window size,
%                     filling missing values from X with -Inf. This is
%                     equivalent to padding X with -Inf at the endpoints.
%       'discard'   - compute the maximum only when the window is filled
%                     with elements of X, discarding partial endpoint
%                     calculations and their corresponding elements in Y.
%                     This truncates the output; for a vector X and window
%                     length K, Y has length LENGTH(X)-K+1.
%                     
%   When ENDPT is a scalar numeric or logical value, the missing elements
%   of X inside the window are replaced with that value and Y remains the
%   same size as X.
%
%   Example: Compute a 5-point centered moving maximum.
%       t = 1:10;
%       x = [4 3 1 -2 -4 -1 0 2 -2 3];
%       yc = movmax(x,5);
%       plot(t,x,t,yc);
%
%   Example: Compute a 5-point trailing moving maximum.
%       t = 1:10;
%       x = [4 3 1 -2 -4 -1 0 2 -2 3];
%       yt = movmax(x,[4 0]);
%       plot(t,x,t,yt);
%
%   Example: Compute a 5-point centered moving maximum, padding the ends of
%   the input with -Inf.
%       t = 1:10;
%       x = [4 3 1 -2 -4 -1 0 2 -2 3];
%       yp = movmax(x,5,'Endpoints','fill');
%       plot(t,x,t,yp);
%
%   Example: Compute a 5-point trailing moving maximum, ignoring the first
%   4 window shifts that do not contain 5 input elements.
%       x = [4 3 1 -2 -4 -1 0 2 -2 3];
%       yd = movmax(x,[4 0],'Endpoints','discard');
%
%   See also MAX, CUMMAX, MOVMIN, MOVSUM
%   

% Copyright 2015 The MathWorks, Inc.
% Built-in function.