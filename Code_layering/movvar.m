%MOVVAR   Moving variance value.
%   Y = MOVVAR(X,K) for a vector X and positive integer scalar K computes a
%   centered moving variance by sliding a window of length K along X. Each
%   element of Y is the local variance of the corresponding values of X
%   inside the window, with Y the same size as X. When K is even, the
%   window is centered about the current and previous elements of X. The
%   sliding window is truncated at the endpoints where there are fewer than
%   K elements from X to fill the window.
%   
%   For N-D arrays, MOVVAR operates along the first array dimension whose
%   size does not equal 1.
%
%   By default, MOVVAR normalizes by K-1 if K>1. If X consists of
%   independent, identically distributed samples, then MOVVAR is an
%   unbiased estimator of the variance of the population of each window.
%
%   Y = MOVVAR(X,[NB NF]) for a vector X and nonnegative integers NB and NF
%   computes a moving variance along the length of X, returning the local
%   variance of the previous NB elements, the current element, and the next
%   NF elements of X.
%
%   MOVVAR(X,K,NRM) specifies the normalization factor for the variance and
%   can be one of the following:
%
%       0   - (default) normalizes by K-1 for K>1 and by K when K=1.
%       1   - normalizes by K and produces the second moment of the window
%             about its mean.
%
%   Y = MOVVAR(X,K,NRM,DIM) or Y = MOVVAR(X,[NB NF],NRM,DIM) operates along
%   dimension DIM of X. When specifying DIM, you must specify NRM.
%
%   MOVVAR(...,MISSING) specifies how NaN (Not-a-Number) values are treated
%   and can be one of the following:
%
%       'includenan'   - (default) the variance of any window containing
%                        NaN values is also NaN.
%       'omitnan'      - the variance of any window containing NaN values
%                        is the variance of all its non-NaN elements. If
%                        all elements are NaN, the result is NaN.
%
%   MOVVAR(...,'Endpoints',ENDPT) controls how the variance is calculated
%   at the endpoints of X, where there are not enough elements to fill the
%   window. ENDPT can be either a scalar numeric or logical value or one of
%   the following:
%
%       'shrink'    - (default) compute the variance over the number of
%                     elements of X that are inside the window, effectively
%                     reducing the window size to fit X at the endpoints.
%       'fill'      - compute the variance over the full window size,
%                     filling missing values from X with NaN. This is
%                     equivalent to padding X with NaN at the endpoints.
%       'discard'   - compute the variance only when the window is filled
%                     with elements of X, discarding partial endpoint
%                     calculations and their corresponding elements in Y.
%                     This truncates the output; for a vector X and window
%                     length K, Y has length LENGTH(X)-K+1.
%                     
%   When ENDPT is a scalar numeric or logical value, the missing elements
%   of X inside the window are replaced with that value and Y remains the
%   same size as X.
%
%   Example: Compute a 5-point centered moving variance.
%       t = 1:10;
%       x = [4 8 6 -1 -2 -3 -1 3 4 5];
%       yc = movvar(x,5);
%       plot(t,x,t,yc);
%
%   Example: Compute a 5-point trailing moving variance.
%       t = 1:10;
%       x = [4 8 6 -1 -2 -3 -1 3 4 5];
%       yt = movvar(x,[4 0]);
%       plot(t,x,t,yt);
%
%   Example: Compute a 5-point centered moving variance, padding the ends
%   of the input with NaN.
%       t = 1:10;
%       x = [4 8 6 -1 -2 -3 -1 3 4 5];
%       yp = movvar(x,5,'Endpoints','fill');
%       plot(t,x,t,yp);
%
%   Example: Compute a 5-point trailing moving variance, ignoring the first
%   4 window shifts that do not contain 5 input elements.
%       x = [4 8 6 -1 -2 -3 -1 3 4 5];
%       yd = movvar(x,[4 0],'Endpoints','discard');
%
%   See also VAR, MOVMEAN, MOVMEDIAN, MOVSTD
%   

% Copyright 2015 The MathWorks, Inc.
% Built-in function.