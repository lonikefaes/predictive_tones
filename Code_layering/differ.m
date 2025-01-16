function f=differ(x,t,m) 
% DIFFER calculates (s-1)-th divided difference of max((t-x)^(m-1),0) where s is the
% number of knots in sequence t (i.e., s=length(t))
%	
% Usage: f=differ(x,t,m);
%
% Inputs:
%   x    dicrete variable
%   t    knot sequence
%   m    order
%
% Output:
%   f    (s-1)-th divided difference of max((t-x)^(m-1),0)
%
% Remark: Sequence t must be ordered. 

% See http://www.ncrg.aston.ac.uk/Projects/BiOrthog/ for more details

s=length(t);pom=0;
if isempty(x)
    f=[];
    return;
else
    if x(2)<x(1)
        x=x(end:-1:1);pom=1; %flip x
    end
end

if s==1
    f=green_fcn(t-x,m-1);
elseif (t(1)==t(s))
    f=(1/factorial(s-1))*prod(m-s+1:m-1)*green_fcn(t(1)-x,m-s);
else
    f=(differ(x,t(2:s),m)-differ(x,t(1:s-1),m))/(t(s)-t(1));
end;

if pom==1
    f=f(end:-1:1); %flip function values back
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
