 % Decide if point P is in triangles E indexing nodes V in 3D
 function [u, v] = point_in_triangle_3d(p,V,F)

 % vectors
 v0 = V(F(:,3),:) - V(F(:,1),:);
 v1 = V(F(:,2),:) - V(F(:,1),:);
 v2 = bsxfun(@minus,p, V(F(:,1),:));
 
 % dot products
 dot00 = dot(v0, v0, 2);
 dot01 = dot(v0, v1, 2);
 dot02 = dot(v0, v2, 2);
 dot11 = dot(v1, v1, 2);
 dot12 = dot(v1, v2, 2);

 % barycentric coordinates
 invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
 u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
 v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;
 
 
