function y = BOLD_voxels2layers(Data,depth_map,K)
% INPUT:
% Data (Voxels x Time-points)
% K - number of depths
% depth_map - voxel map of normalized corical depths (equivolume or equidistant)

% OUTPUT:
% y - (time-points x depths)

T         = size(Data,2);
layer_ax = linspace(0,1,2*(K)+1);
layer_ax = layer_ax(2:2:end);

[dist ind] = sort(depth_map);

ya = flipud(Data(ind,:));



H = [];
layer_axr = linspace(0,1,2*(K)+1);
layer_axr = layer_axr(2:2:end);

layer_axr_ex = [layer_axr(1)-(layer_axr(2)-layer_axr(1)),...
                layer_axr,layer_axr(end)+(layer_axr(2)-layer_axr(1))];

% create bell-shaped basis function to do voxel averaging coresponding to
% particular depth
for i = 1:length(layer_axr_ex),
   H(:,i) = (exp(-(length(layer_axr_ex).*(dist-layer_axr_ex(i))).^2)).^3;
end

H = H./repmat(sum(H,2),1,length(layer_axr_ex));
H = H(:,2:end-1);

Hs  = H./repmat(sum(H,1),size(H,1),1);
%figure; for i = 1:size(Hs,2), plot(Hs(:,i)); hold on; end
y = (Hs'*ya)';
