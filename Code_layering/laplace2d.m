function [U_out] = laplace2d(U)
[n_x, n_y] = size( U );
U_out = U;

% Step 2
u_to_w = zeros( n_x, n_y );
w_to_u = zeros( 2, n_x * n_y );
m = 0;

for ix = 1:n_x
    for iy = 1:n_y
        if U(ix, iy) == -Inf
            m = m + 1;
            u_to_w(ix, iy) = m;
            w_to_u(:, m) = [ix, iy]';
        end
    end
end

% Create the sparse system of linear equations
M = spalloc( m, m, 5*m );
b = zeros( m, 1 );
for k = 1:m
    % Get the coordinates of the kth point
    xy = w_to_u(:,k);
    % 8 points around
    XY = [xy(1)-1,xy(2);
        xy(1)+1,xy(2),;
        xy(1),xy(2)-1;
        xy(1),xy(2)+1;];
    
    M(k,k) = -4;
    
    for i =1:4
        mi    = u_to_w(XY(i,1),XY(i,2));
        kind  = U(XY(i,1),XY(i,2));
        adj(i) = kind;
        if isnan(kind)
            M(k,k) = M(k,k) + 1;
        end
        if kind==-Inf
            if mi>0
                M(k,mi) = 1;
            end
        end
        
    end
    adj(adj==-Inf)  = 0;
    adj(isnan(adj)) = 0;
    b(k) = -sum(adj);
    
    
end

w = M \ b;

for k = 1:m
    % Copy the value from w into U_out
    U_out(w_to_u(1,k),w_to_u(2,k)) = w(k);
end
