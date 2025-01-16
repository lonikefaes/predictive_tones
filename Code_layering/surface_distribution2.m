function [EVsurf, EDsurf, FV_EV, FV_ED, Vol,PV_EV] = surface_distribution(FVin,fiber,layers,EV_vox)


vertices = FVin.vertices;
faces    = FVin.faces;
% 
% 
% edges   = meshEdges(faces);
% Ne      = vertex_neighbours(FVin);
xi      = linspace(0,1,layers);
% xi(1)   = -1/(layers-1)*0.015; 
% xi(end) =  1+1/(layers-1)*0.065; 

% r2      = cell(size(vertices,1),1); 
% [r2{:}] = deal(ones(99,20)*NaN);

%N  = patchnormals(FVin);

%   
%
%   figure, 
%     patch(FVin,'facecolor',[0 0 1]); view(3); hold on
%    [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FVin);
%     p1=FVin.vertices; p2=FVin.vertices+10*diag(Lambda1)*Dir1;       
%     plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
%     axis equal; view(3); 
% 
%     p1=FVin.vertices; p2=FVin.vertices+10*diag(Lambda1)*Dir2;       
%     plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','r-');
%     axis equal; view(3); hold off
close all;
figure; hold on
for l =1:size(fiber{1},1)
    FVx = FVin;
    for i = 1:length(fiber)
        FVx.vertices(i,:) = fiber{i}(l,:);
    end
    [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FVx);
    p1=FVx.vertices(1:10:200,:); p2=FVx.vertices(1:10:200,:)+10*diag(sqrt(abs(Lambda1(1:10:200))))*Dir1(1:10:200,:);
    plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
    axis equal; view(3); hold on
    p1=FVx.vertices(1:10:200,:); p2=FVx.vertices(1:10:200,:)+10*diag(sqrt(abs(Lambda2(1:10:200))))*Dir2(1:10:200,:);
    plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','r-');
    axis equal; view(3); drawnow;
    AA(:,l) = sqrt(abs(Lambda1)).*sqrt(abs(Lambda2));
end
%Cgaussian = Cgaussian + 1;
%V = cell(size(vertices,1),1); 
for j = 1:length(fiber)
    h(j,:) = sqrt(sum(diff(fiber{j},[],1).^2,2));
end
for i = 1:size(fiber{1}),
    
  
    FVx = FVin;
    for j = 1:length(fiber)
        FVx.vertices(j,:) = fiber{j}(i,:);

    end
     [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(FVx,1);
    % p1=FVx.vertices(1000,:); p2=FVx.vertices(1000,:)+20*diag(sqrt(abs(Lambda1(1000))))*Dir1(1000,:);
    % plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','g-');
    % axis equal; view(3); hold on
    % p1=FVx.vertices(1000,:); p2=FVx.vertices(1000,:)+20*diag(sqrt(abs(Lambda2(1000))))*Dir2(1000,:);
    % plot3([p1(:,1) p2(:,1)]',[p1(:,2) p2(:,2)]',[p1(:,3) p2(:,3)]','r-');
    % axis equal; view(3); drawnow;
    A(:,i) = sqrt(abs(Lambda1)).*sqrt(abs(Lambda2));
   
  %  test  = interp1(V{i},V{i},xi,'pchip','extrap');
end;

  %  V{i}  = pi.*d2{i}.*(r1{i}(1:end-1).^2 + r1{i}(1:end-1).*r1{i}(2:end) + r1{i}(2:end).^2)/3;
    V     = h.*pi.*(A(:,2:end) + A(:,1:end-1) + sqrt(A(:,2:end).*A(:,1:end-1)));
     V  = cumsum(V,2); 
     V  = V - repmat(min(V,[],2),1,size(V,2));
     V  = V./repmat(max(V,[],2),1,size(V,2));
% 
     D  = cumsum(h,2); 
     D  = D - repmat(min(D,[],2),1,size(D,2));
     D  = D./repmat(max(D,[],2),1,size(D,2));
%     
 for j = 1:length(fiber)   
    EVsurf{j}  = interp1(V(j,:)',fiber{j}(2:end,:),xi,'pchip','extrap');
    EDsurf{j}  = interp1(D(j,:)',fiber{j}(2:end,:),xi,'pchip','extrap');
    Vol{j}     = xi(:);
 end

for j = 1:layers,
    for i = 1:length(fiber),
        FV_EV{j,1}.vertices(i,:) = EVsurf{i}(j,:);
        FV_ED{j,1}.vertices(i,:) = EDsurf{i}(j,:);
    end;
    FV_EV{j,1}.faces = faces;
    FV_ED{j,1}.faces = faces;
end;

if ~isempty(EV_vox)

FV_EV = surface_smooth(FV_EV,1);
% iterate
figure,
xii = xi; 
Iter = 5;
cm = jet(Iter);
for it = 1:Iter
    disp(it),
    PV_EV = surface_voxel_partial_volume_simple(FV_EV,EV_vox,Vol{1},0);
    vol_EV = [];
    for i = 1:layers-1
        vol_EV(i) = sum(vec(PV_EV{i}(:,:,6:size(EV_vox,3)-5)));
    end;
    plot([vol_EV'],'color',cm(it,:)); hold on; drawnow;
    
    m = median(vol_EV);
    
    Norm = vol_EV;
    Norm = Norm - m;
    %if it == 1;
       maxN = max(abs(Norm)); 
   % end
    Norm = Norm./maxN;
    xii = [0,cumsum(diff(xii)-Norm*(0.001))];
    xii = xii./max(xii);
    xii = (1.01+0.01)*xii - 0.01
    for i = 1:length(fiber),
        EVsurf{i}  = interp1(V{i},fiber{i}(2:end,:),xii,'pchip','extrap');
    end
    
    for j = 1:layers,
        for i = 1:length(fiber),
            FV_EV{j,1}.vertices(i,:) = EVsurf{i}(j,:);
        end;
        FV_EV{j,1}.faces = faces;
    end;
end

else
    PV_EV = [];
end