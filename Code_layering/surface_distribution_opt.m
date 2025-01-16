function [EVsurf, EDsurf, FV_EV, FV_ED, FV_EVmid, FV_EDmid, Vol] = surface_distribution_opt(FVin,fiber,layers,Iter)

if nargin<4
   Iter = 0; 
end

vertices = FVin.vertices;
faces    = FVin.faces;

TR           = triangulation(faces,vertices);
[~,FBpoints] = freeBoundary(TR);
FBtri = find(ismember(vertices,FBpoints,'rows')==1);

Ne      = vertex_neighbours(FVin);
xi      = linspace(0,1,layers);
% xi(1)   = -1/(layers-1)*0.005; 
% xi(end) =  1+1/(layers-1)*0.03; 

r2      = cell(size(vertices,1),1); 
[r2{:}] = deal(ones(99,20)*NaN);
%FVx  = FVin;
% for j = 1:size(fiber{1},1)
%     for i = 1:length(fiber),
%         FVx.vertices(i,:) = fiber{i}(j,:);
%         P{j}(i,:) =  fiber{i}(j,:);
%     end
%     N{j} = patchnormals(FVx);
%     plane{j} = createPlane(P{j},N{j});
% end
for i = 1:length(fiber),
%     for j = 1:size(fiber{1},1)
%         NN{i}(j,:) = N{j}(i,:);
%         PP{i}(j,:) = plane{j}(i,:);
%     end
    nn2 = fiber{i}(1:end-1,:) - fiber{i}(2:end,:);
    nn2 = nn2./repmat(sqrt(sum(nn2.^2,2)),1,3);
    NN2{i} = interp1([1:size(nn2)]',nn2,[0:size(nn2)]','linear','extrap');
    PP2{i} = createPlane(fiber{i},NN2{i});
    P(i,:) = fiber{i}(round(size(fiber{i},1)/2),:);
end
%   

V = cell(size(vertices,1),1); 

for i = 1:length(fiber),
    
    N1   = unique(Ne{i}); % fisrt ring
    N12  = unique([Ne{Ne{i}}]); % first and second ring
    N123 = unique([Ne{[Ne{Ne{i}}]}]);
    N12(N12==i) = [];
    N1(N1==i) = [];
    N123(N123==i) = [];
   % N2  = N12(~ismember(N12,N1)); % second ring;
    if ~logical(sum(ismember(FBtri,i)))
        for j = 1:length(N12)
            if ~logical(sum(ismember(FBtri,N12(j))));
                ProjP = projPointOnPlane(fiber{N12(j)},PP2{i});
                r2{i}(:,j) = sqrt(sum((fiber{i} - ProjP).^2,2));
            end
        end
    else
           FBneighbors{i} = N123(~ismember(N123',FBtri));
    end



  %  r2{i} = (1-0.3)*nanmean(r1{i},2) + 0.3*nanmean(r3{i},2);
    d2{i} = sqrt(sum(diff(fiber{i},[],1).^2,2));


    
    if ~logical(sum(ismember(FBtri,i)))
        r2{i} = nanmean(r2{i},2);
        V0{i}  = pi.*d2{i}.*(r2{i}(1:end-1).^2 + r2{i}(1:end-1).*r2{i}(2:end) + r2{i}(2:end).^2)/3;
        As    = pi*r2{i}(2:end).^2;
        Al    = pi*r2{i}(1:end-1).^2;
        V{i}  = d2{i}./4.*(As + (As.^2.*Al).^(1/3) + (As.*Al.^2).^(1/3) + Al);
        V{i}  = cumsum(V{i});
        V{i}  = V{i} - min(V{i});
        V{i}  = V{i}./max(V{i});
        
        V0{i}  = cumsum(V0{i});
        V0{i}  = V0{i} - min(V0{i});
        V0{i}  = V0{i}./max(V0{i});
     %   figure(1), plot([V0{i},V{i}]);
        EVsurf{i}  = interp1(V{i},fiber{i}(2:end,:),xi,'pchip','extrap');
    end
    
    D{i}  = cumsum(d2{i});
    D{i}  = D{i} - min(D{i});
    D{i}  = D{i}./max(D{i});
    EDsurf{i}  = interp1(D{i},fiber{i}(2:end,:),xi,'pchip','extrap');
 
    Vol{i}     = xi(:);
  %  test  = interp1(V{i},V{i},xi,'pchip','extrap');
end;

for i = 1:size(FBtri,1)
    Vx = [];
    for j = 1:length(FBneighbors{FBtri(i)})
        if j==1
            Vx = V{FBneighbors{FBtri(i)}(j)};
        else
            Vx = Vx+V{FBneighbors{FBtri(i)}(j)};
        end
    end
    Vx = Vx/j;
    V{FBtri(i)} = Vx;
    EVsurf{FBtri(i,1)}  = interp1(Vx,fiber{FBtri(i)}(2:end,:),xi,'pchip','extrap');
end

for j = 1:layers,
    
    for i = 1:length(fiber),
        FV_EV{j,1}.vertices(i,:) = EVsurf{i}(j,:);
        FV_ED{j,1}.vertices(i,:) = EDsurf{i}(j,:);
        if j<layers
            FV_EVmid{j,1}.vertices(i,:) = (EVsurf{i}(j,:) + EVsurf{i}(j+1,:))/2;
            FV_EDmid{j,1}.vertices(i,:) = (EDsurf{i}(j,:) + EDsurf{i}(j+1,:))/2;
        end
    end;

    FV_EV{j,1}.faces = faces;
    FV_ED{j,1}.faces = faces;
     if j<layers
         FV_EVmid{j,1}.faces = faces;
         FV_EDmid{j,1}.faces = faces;
     end
end




if Iter>1
xii = xi; 
FV_EV = surface_smooth(FV_EV,1);
cm = jet(Iter);
for it = 1:Iter
    disp(it),
    PV_EV    = partial_volume(FV_EV);

    plot([PV_EV'],'color',cm(it,:)); hold on; drawnow;
    
    m = median(PV_EV);
    
    Norm = PV_EV;
    Norm = Norm - m;
    if it == 1;
       maxN = std(abs(Norm)); 
    end
    Norm = Norm./maxN;
    xii = [0,cumsum(diff(xii)-Norm*(0.001))];
    xii = xii./max(xii);
    xii = (1.01+0.01)*xii - 0.01;
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PV = partial_volume(FVs),

faces  = FVs{1,1}.faces;
layers = size(FVs,1);
V      = zeros(size(faces,1),layers-1);

for i = 1:size(faces,1)
    
    ind = faces(i,:);
    for j = 1:layers;
        Pts     = FVs{j}.vertices(ind,:);
        A(j,1)  = triangleArea3d(Pts);
        P1(j,:) = Pts(1,:);
        P2(j,:) = Pts(2,:);
        P3(j,:) = Pts(3,:);
    end;
    
    D1 = sqrt(sum(diff(P1).^2,2));
    D2 = sqrt(sum(diff(P2).^2,2));
    D3 = sqrt(sum(diff(P3).^2,2));
    
    V1 = D1.*(A(1:end-1) + sqrt(A(1:end-1).*A(2:end)) + A(2:end))./3;
    V2 = D2.*(A(1:end-1) + sqrt(A(1:end-1).*A(2:end)) + A(2:end))./3;
    V3 = D3.*(A(1:end-1) + sqrt(A(1:end-1).*A(2:end)) + A(2:end))./3;
    
    V(i,:) = (V1+V2+V3)./3;
end;

PV = sum(V,1);


% 
% if ~isempty(EV_vox)
% 
% FV_EV = surface_smooth(FV_EV,1);
% % iterate
% figure,
% xii = xi; 
% Iter = 5;
% cm = jet(Iter);
% for it = 1:Iter
%     disp(it),
%     PV_EV = surface_voxel_partial_volume_simple(FV_EV,EV_vox,Vol{1},0);
%     vol_EV = [];
%     for i = 1:layers-1
%         vol_EV(i) = sum(vec(PV_EV{i}(:,:,6:size(EV_vox,3)-5)));
%     end;
%     plot([vol_EV'],'color',cm(it,:)); hold on; drawnow;
%     
%     m = median(vol_EV);
%     
%     Norm = vol_EV;
%     Norm = Norm - m;
%     %if it == 1;
%        maxN = max(abs(Norm)); 
%    % end
%     Norm = Norm./maxN;
%     xii = [0,cumsum(diff(xii)-Norm*(0.001))];
%     xii = xii./max(xii);
%     xii = (1.01+0.01)*xii - 0.01;
%     for i = 1:length(fiber),
%         EVsurf{i}  = interp1(V{i},fiber{i}(6:end-5,:),xii,'linear','extrap');
%     end
%     
%     for j = 1:layers,
%         for i = 1:length(fiber),
%             FV_EV{j,1}.vertices(i,:) = EVsurf{i}(j,:);
%         end;
%         FV_EV{j,1}.faces = faces;
%     end;
% end
% 
% else
%     PV_EV = [];
% end