function FVout = surface_smooth(FVin,iter),


if iscell(FVin)
     FVout = cell(size(FVin,1),1);
%     [FVout{:}] = deal(struct('vertices',FVin{1}.vertices*0,'faces',FVin{1}.faces));
    parfor i = 1:size(FVin,1)

        [FVout{i,1}] = smoothMesh(FVin{i}.vertices, FVin{i}.faces,iter,1);
    end
else
    [FVout] = smoothMesh(FVin.vertices, FVin.faces,iter,1);
end