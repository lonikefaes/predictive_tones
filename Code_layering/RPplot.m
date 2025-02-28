function [RP,DD,X] = RPplot(P,m,t,r,index)

%  Recurrence plot (RP) with fixed threshold value 

%  Input:   P: time series;
%           m: the embedding dimension
%           t: the delay time
%           r: the threshold value 
%           index: the value 0 or 1 (0 for Euclidean distance,1 for maximum distance)
% Output: 
%           RP:  Recurrence plot Matrix
%           DD:  Distance Matrix          

%  referrence: 
%  X Li, G Ouyang, X Yao, X Guan, Dynamical characteristics of pre-epileptic seizures in rats with recurrence quantification analysis
%  Physics Letters A 333 (1), 164-171

%%  If you need these codes that implement critical functions with (fast) C code, please visit my website:
%%  http://www.escience.cn/people/gxouyang/Tools.html

%  revise time: May 5 2014, Ouyang,Gaoxiang
%  Email: ouyang@bnu.edu.cn
%  
%  Example:
% P = randn(1,100);
% [RP,DD] = RPplot(P,3,1,.5,0);
% subplot(121);imagesc(RP)
% P = sin(0.1:0.1:10);
% [RP,DD] = RPplot(P,3,1,.5,0);
% subplot(122);imagesc(RP)

r = r*std(P);

if nargin < 4, error('Not enough input arguments'),end
N=length(P);

X=zeros(N-(m-1)*t,m);
for k=1:(N-(m-1)*t)
    PP=[];
    for i=1:m
        PP=[PP P(1,k-t+i*t)];
    end
    X(k,:)=PP;
end
N1=length(X);

DD=zeros(N1,N1);
if index == 0;
    r = r*r;
    for k1=1:N1 
       for k2=k1+1:N1
           DD(k1,k2) = sum( (X(k1,:) - X(k2,:)).^2 );
       end
     end
    DD=DD+DD';
end

if index == 1
    for k1=1:N1       
        for k2=k1+1:N1
            DD(k1,k2) = max(abs(X(k1,:)-X(k2,:)));
        end
    end
    DD=DD+DD';
end

RP=zeros(N1,N1);
for i=1:N1       
    for j=i+1:N1
        if DD(i,j) <= r
            RP(i,j)=1;
            RP(j,i)=1;
        end
    end
end
RP = RP +eye(N1);


