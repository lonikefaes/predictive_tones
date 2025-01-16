function segement_cluster(Smoothed,SegIndex)

Position = [1:size(Smoothed,1)];
len = 20;
%for iloop = 1:length(GM05296_Data)
    seg_num = numel(SegIndex) - 1;
    if seg_num > 1
        
        segidx       = SegIndex;
        segidx_emadj = SegIndex;
        
        for j = 2:seg_num
            ileft = min(segidx(j) - len, segidx(j));
            iright = max(segidx(j) + len, segidx(j));
            gmx = Position(ileft:iright);
            gmy = Smoothed(ileft:iright);
            
            % Select initial guess for the of cluster index for each point.
            gmpart = (gmy > (min(gmy) + range(gmy)/2)) + 1;
            
            % Create a Gaussian mixture model object
            gm = gmdistribution.fit(gmy, 2, 'start', gmpart);
            gmid = cluster(gm,gmy);
            
            segidx_emadj(j) = find(abs(diff(gmid))==1) + ileft;
            
        end
        
        % Remove repeat indices
        zeroidx  = [diff(segidx_emadj) == 0; 0];
        SegIndex = segidx_emadj(~zeroidx);
    end
