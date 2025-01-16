function [Ts] = glm_voxel(Y1,Y2)

[x y z T] = size(Y1);



Y1 = reshape(Y1,[x*y*z,T])';
Y2 = reshape(Y2,[x*y*z,T])';
baseline = [1:4];
positive = [7:11];
negative = [14:16];

Y = [Y1(baseline,:);Y1(positive,:);Y1(negative,:);Y2(baseline,:);Y2(positive,:);Y2(negative,:);];

reg1 = zeros(size(Y,1)/2,1);
reg2 = zeros(size(Y,1)/2,1);
reg1([length(baseline)+1:length(baseline)+length(positive)]) = 1;
reg2([length(baseline)+length(positive):size(Y,1)/2]) = 1;
X = [kron(eye(2),[reg1,reg2]),kron(eye(2),ones(size(Y,1)/2,1))];
    
con{1} = [1 0 0 0 0 0]'; % positive
con{2}=  [0 1 0 0 0 0]';
con{3}  = [0 0 1 0 0 0]';
con{4}  = [0 0 0 1 0 0]';
con{5} = [1 0 -1 0 0 0]';
con{6}  = [0 1 0 -1 0 0]';
ts = zeros(size(Y,2),length(con));
for ii = 1:size(Y,2),
    b  = X\Y(:,ii);
    e  = Y(:,ii) - X*b;
    s2 = (e'*e)/(size(X,1)-size(X,2));
    for c = 1:length(con),
        ts(ii,c) =  (con{c}'*b)./sqrt(s2*con{c}'*inv(X'*X)*con{c});
    end
end

Ts = reshape(ts,[x y z c]);
