function M = region_split(Mask_morph)

[nx,ny,nz] = size(Mask_morph); 

[Label, L]  = bwlabeln(Mask_morph==1);

no_regions  = L;

if no_regions>1
    
    for i = 1:no_regions
        BW{i}         = bwdist(Label==i);
    end
    M = cell(no_regions,1);
    [M{:}] = deal(zeros(nx,ny,nz));
    for i = 1:no_regions
        comp    = [1:no_regions];
        comp(i) = [];
        M{i}    = BW{i};
        if no_regions>1
            for k = comp,
                M{i} = M{i}< BW{k};
            end
        end
        % M{i} = M{i}*0;
        %  M{i}(M{i}>0) = 0;
        M{i}(M{i}<0) = 0;
        if  no_regions==1
            M{i} = ~M{i};
            M{i} = imdilate(M{i},nhood);
        end
        M{i}         = logical(M{i});
    end
else
    % if only one region
    M{1} = ones(nx,ny,nz);
end
