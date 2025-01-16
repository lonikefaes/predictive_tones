function POLY = surface2polylines(LS,nz,M,interp_factor)
no_regions = length(M);

for r = 1:no_regions
    for k = 1:nz
        for l = 1:length(LS)
            P = intersectPlaneSurf(LS{l,r},[0 0 (k*interp_factor)-interp_factor/2],[0 0 1]);
            POLY{k,l,r} = P{1}';
        end
        disp(k);
    end
    
    
end