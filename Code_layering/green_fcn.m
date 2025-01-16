function f=green_fcn(x,m)
% GREEN calculates (x_+)^m = x^m for x>=0 (it is 0 for x<0). This functions is known as 'truncated powers'.
%
% Usage: f=green(x,m);
%
% Inputs:
%   x    discrete variable
%   m    order
%
% Output:
%   f    (x_+)^m (see definition of this above)
%
% References:
%    L.L. Schumaker, Spline Functions: Basic Theory, New York-Wiley, 1981.

% See http://www.ncrg.aston.ac.uk/Projects/BiOrthog/ for more details

if m>0
    f=max(x,0).^m;
elseif m==0
    x(x>0)=1;x(x<=0)=0;
    f=x;
else f=0; %too many derivatives in differ.m
end;


%Copyright (C) 2006 Miroslav ANDRLE and Laura REBOLLO-NEIRA
%
%This program is free software; you can redistribute it and/or modify it under the terms 
%of the GNU General Public License as published by the Free Software Foundation; either 
%version 2 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
%without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program;
%if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
%Boston, MA  02110-1301, USA.
