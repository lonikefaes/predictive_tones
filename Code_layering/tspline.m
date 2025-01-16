function D=tspline(x,t,m),
% TSPLINE generates B-spline basis of order m corresponding to knot sequence t
% 
% Usage:  D=tspline(x,t,m);
%           
% Inputs:
%   x     variable range (discrete points)
%   t     knot sequence (could be multiple knots)
%   m     order of splines (m=1 means a piecewise constant functions)  
% Output:
%   D     B-spline function basis corresponding to knot sequence t
%  
% Remark: The sequence t is sorted before calculation.
%  
% References: 
%    L.L. Schumaker, Spline Functions: Basic Theory, New York, Wiley, 1981.  

% See http://www.ncrg.aston.ac.uk/Projects/BiOrthog/ for more details
 
l=length(t);
s=length(x);
% sorting the knot sequence
t=sort(t); 

if m<1
  error('Order of splines, m, has to be >=1');
else  
  if l<m+1 
    error('Too short knot sequence, the minimal length is m+1 points');
  end;
  D=zeros(s,l-m);
  for i=1:l-m
    D(:,i)=[(t(m+i)-t(i))*differ(x,t(i:i+m),m)'];
  end  
end


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
