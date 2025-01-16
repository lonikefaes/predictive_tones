function out = isPointInTriangle(P,V,F),

 nPts = size(P,1);
 nTri = size(F,1);
 
 u = zeros(nPts, nTri);
 v = zeros(nPts, nTri);
 for i = 1:nPts
     [u(i,:), v(i,:)] = point_in_triangle_3d(P(i,:),V,F);
 end
 
out = u >= -eps & v >= -eps & (u+v-eps) <= 1; 