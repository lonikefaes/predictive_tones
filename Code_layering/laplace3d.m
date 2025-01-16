function [U_out] = laplace3d(U)

U_out               = U;
greyMatter          = find(U_out == -Inf);
N                   = length(greyMatter);

u_to_w              = zeros(size(U_out));
u_to_w(greyMatter)  = 1:N;
[x, y, z]           = ind2sub(size(U_out), greyMatter);
w_to_u              = [x, y, z];

XYZ          = repmat(w_to_u, [1, 1, 6]);
XYZ(:, 1, 1) = XYZ(:, 1, 1) - 1;
XYZ(:, 1, 2) = XYZ(:, 1, 2) + 1;
XYZ(:, 2, 3) = XYZ(:, 2, 3) - 1;
XYZ(:, 2, 4) = XYZ(:, 2, 4) + 1;
XYZ(:, 3, 5) = XYZ(:, 3, 5) - 1;
XYZ(:, 3, 6) = XYZ(:, 3, 6) + 1;

mi      = squeeze(u_to_w(sub2ind(size(U_out), XYZ(:, 1, :), XYZ(:, 2, :), XYZ(:, 3, :))));
kind    = squeeze(U_out(sub2ind(size(U_out), XYZ(:, 1, :), XYZ(:, 2, :), XYZ(:, 3, :))));

indices     = kind == -Inf & mi > 0;
k           = mod(find(indices), N);
k(k == 0)   = N;

i = [(1:N)'; k(:)];
j = [(1:N)'; mi(indices)];
v = [-6 * ones(N, 1) + sum(isnan(kind), 2); ones(length(k), 1)];
M = sparse(i, j, v, N, N);

kind(kind == -Inf)  = 0;
kind(isnan(kind))   = 0;
b                   = -sum(kind, 2);  

U_out(sub2ind(size(U_out), w_to_u(:, 1), w_to_u(:, 2), w_to_u(:, 3))) = M \ b;



return
[n_x, n_y, n_z] = size( U );
    U_out_dx = U*0;
    N = length(find(U==-Inf));
    % Step 2
    u_to_w = zeros( n_x, n_y, n_z );
    w_to_u = zeros( 3,N);
    m = 0;
    
    for iz = 1:n_z
        for iy = 1:n_y
            for ix = 1:n_x
                if U(ix, iy, iz) == -Inf
                    m = m + 1;
                    u_to_w(ix, iy, iz) = m;
                    w_to_u(:, m) = [ix, iy, iz]';
                end
            end
        end
    end

    % Create the sparse system of linear equations
    M = spalloc( m, m, 7*m );
    b = zeros( m, 1 );


    adj = [];
    for k = 1:m
        % Get the coordinates of the kth point
        xyz = w_to_u(:,k);
        % 6 points around
        XYZ = [xyz(1)-1,xyz(2),xyz(3);
               xyz(1)+1,xyz(2),xyz(3);];
             %  xyz(1),xyz(2)-1,xyz(3);
              % xyz(1),xyz(2)+1,xyz(3);
              % xyz(1),xyz(2),xyz(3)-1;
              % xyz(1),xyz(2),xyz(3)+1;];
          
        M(k,k) = -2;
        F(k,1)= U_out(xyz(1),xyz(2),xyz(3));%./max(U_out(:));
        for i =1:2
            mi    = u_to_w(XYZ(i,1),XYZ(i,2),XYZ(i,3));
            kind  = U(XYZ(i,1),XYZ(i,2),XYZ(i,3));
            
            adj(i) = 0;
            if isnan(kind)
                M(k,k) = M(k,k) + 1;
            
            elseif kind==-Inf
                if mi>0
                    M(k,mi) = 1;
                end
            else
                count(i) = 1
            end
           
          
        end
        adj(adj==-Inf)  = 0;
        adj(isnan(adj)) = 0;
        b(k) = -sum(adj);
 
      
    end

    w = M \(F+b);

    for k = 1:m
        % Copy the value from w into U_out
        U_out_dx(w_to_u(1,k),w_to_u(2,k),w_to_u(3,k)) = w(k);
    end
end