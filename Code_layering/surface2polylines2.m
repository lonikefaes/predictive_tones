function POLY = surface2polylines2(LS,nz,interp_factor)




            P = intersectPlaneSurf(LS,[0 0 (nz*interp_factor)-interp_factor/2],[0 0 1]);
          %  POLY= P{1}';
         [POLY(:,1),POLY(:,2)] = poly2cw(P{1}(1,:),P{1}(2,:));

