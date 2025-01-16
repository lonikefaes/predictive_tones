function [EVsurf, EDsurf, FV_EV, FV_ED, Vol,PV_EV, Thickness,Curvature,Thickness2,CosAngle] = surface_distribution(FVin,fiber,layers,EV_vox)

if nargin<4
   EV_vox = []; 
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
FVx  = FVin;
for j = 1:size(fiber{1},1)
    for i = 1:length(fiber),
        FVx.vertices(i,:) = fiber{i}(j,:);
      %  P{j}(i,:) =  fiber{i}(j,:);
    end
    N{j} = patchnormals(FVx);
   % plane{j} = createPlane(P{j},N{j});
end
for i = 1:length(fiber),
    for j = 1:size(fiber{1},1)
        NN{i}(j,:) = N{j}(i,:);
      %  PP{i}(j,:) = plane{j}(i,:);
    end
    nn2 = fiber{i}(1:end-1,:) - fiber{i}(2:end,:);
    nn2 = nn2./repmat(sqrt(sum(nn2.^2,2)),1,3);
    NN2{i} = interp1([1:size(nn2)]',nn2,[0:size(nn2)]','linear','extrap');
    PP2{i} = createPlane(fiber{i},NN2{i});
    P(i,:) = fiber{i}(round(size(fiber{i},1)/2),:);
end
%   

V = cell(size(vertices,1),1); 
DD = cell(size(vertices,1),1); 
D = cell(size(vertices,1),1); 
%C = cell(size(vertices,1),1); 
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


u  = repmat([0 1 0],size(NN{i},1),1);
v  = [NN{i}];
%Theta{i}= atan2d(norm(cross(u,v)),dot(u',v')');
for j = 1:size(u,1)
Theta2{i}(j,1)= abs(v(j,:)*u(j,:)'./(norm(u(j,:))*norm(v(j,:))));
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
        C{i} = mean(As(1:20))./mean(As(end-20:end));
        %figure(1), plot([mean(As(1:20))./mean(As(end-20:end))],'*');
        EVsurf{i}  = interp1(V{i}(5:end-5,:),fiber{i}(6:end-5,:),xi,'linear','extrap');
        Curvature{i} = interp1(V{i}(5:end-5,:),V{i}(5:end-5)./V{i}(5:end-5).*C{i},xi,'linear','extrap');
       % pause
    end
    
    D{i}  = cumsum(d2{i});
    D{i}  = D{i} - min(D{i});
    D{i}  = D{i}./max(D{i});
    EDsurf{i}  = interp1(D{i}(5:end-5,:),fiber{i}(6:end-5,:),xi,'linear','extrap');
    DD{i}      = cumsum(d2{i});
    Thickness{i} = interp1(D{i}(5:end-5,:),DD{i}(5:end-5),xi,'linear','extrap');
    CosAngle{i} = interp1(xi(1:end-1),Theta2{i},xi,'linear','extrap');
    Thickness2{i} = interp1(D{i}(5:end-5,:),DD{i}(5:end-5)./DD{i}(5:end-5)*max(DD{i}(5:end-5)),xi,'linear','extrap');
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
    EVsurf{FBtri(i,1)}  = interp1(Vx(5:end-5,:),fiber{FBtri(i,1)}(6:end-5,:),xi,'linear','extrap');
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
    xii = (1.01+0.01)*xii - 0.01;
    for i = 1:length(fiber),
        EVsurf{i}  = interp1(V{i},fiber{i}(6:end-5,:),xii,'linear','extrap');
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