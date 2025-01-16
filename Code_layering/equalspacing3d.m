function [Xa Ya Za] = equalspacing3d(xyz0,up,flag,norm,th),

if nargin<4
    norm = [];
end
if nargin<5
   th = 3; 
end
xyz0(isnan(xyz0(:,1)),:) = [];
%check for break
ls = sqrt(sum(diff(xyz0,[],1).^2,2));
breaks = find(abs(ls)>th);
breaks = [0; breaks(:); size(xyz0,1)];
upt = diff(breaks);

if nargin<3
    flag = 0;
end
if flag == 1
    up = up.*upt./breaks(end);
end
Xa = []; Ya = []; Za = [];
for i =1:length(breaks)-1
    
    xyz = xyz0(breaks(i)+1:breaks(i+1),:);


 %   if down~=1
%         X = spm_interp(xyz(:,1),down);
%         Y = spm_interp(xyz(:,2),down);
%         Z = spm_interp(xyz(:,3),down);
 %   else
        X = xyz(:,1);
        Y = xyz(:,2);
        Z = xyz(:,3);
   % end
    
%     if length(X)<length(xyz(:,1))
%         X = [xyz(1,1);X];
%         Y = [xyz(1,2);Y];
%         Z = [xyz(1,3);Z];
%     end
  %  if isempty(norm)
    norm = cumsum([0;sqrt(sum(diff([X,Y,Z],[],1).^2,2))+1e-8]);
   % end
    if ~flag,
        if up~=1
            Xi = interp1(norm,X,[0:1/up:norm(end)]','cubic');
            Yi = interp1(norm,Y,[0:1/up:norm(end)]','cubic');
            Zi = interp1(norm,Z,[0:1/up:norm(end)]','cubic');
        else
            Xi = X;
            Yi = Y;
            Zi = Z;
        end
    else
        if up(i)~=1
            Xi = interp1(norm,X,linspace(0,norm(end),up(i))','cubic');
            Yi = interp1(norm,Y,linspace(0,norm(end),up(i))','cubic');
            Zi = interp1(norm,Z,linspace(0,norm(end),up(i))','cubic');
        else
            Xi = X;
            Yi = Y;
            Zi = Z;
        end
    end
    if length(breaks)>2
        Xa = [Xa; Xi; NaN];   Ya = [Ya; Yi; NaN];   Za = [Za; Zi; NaN];
    else
        Xa = Xi;   Ya = Yi;   Za = Zi;
    end
     
end

if nargout==1
    Xa = [Xa,Ya,Za];
end