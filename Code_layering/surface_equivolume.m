function [EVsurf, EDsurf, FV_EV, FV_ED] = surface_equivolume(FVin,fiber,lam1_in,lam2_in,lam1_out,lam2_out,layers,EV_vox)

if nargin<8
   EV_vox = []; 
end
vertices = FVin.vertices;
faces    = FVin.faces;
% 
% TR           = triangulation(faces,vertices);
% [~,FBpoints] = freeBoundary(TR);
% FBtri = find(ismember(vertices,FBpoints,'rows')==1);
% 
% Ne      = vertex_neighbours(FVin);
% xi      = linspace(0,1,layers);
% xi(1)   = -1/(layers-1)*0.005; 
% xi(end) =  1+1/(layers-1)*0.03; 
% 
% r2      = cell(size(vertices,1),1); 
% [r2{:}] = deal(ones(99,20)*NaN);
%FVx  = FVin;
% for j = 1:size(fiber{1},1)
%     for i = 1:length(fiber),
%         FVx.vertices(i,:) = fiber{i}(j,:);
%         P{j}(i,:) =  fiber{i}(j,:);
%     end
%     N{j} = patchnormals(FVx);
%     plane{j} = createPlane(P{j},N{j});
% end
   

%V = cell(size(vertices,1),1);
  figure(1);
for i = 1:length(fiber),
    
    
    d{i} = sqrt(sum(diff(fiber{i},[],1).^2,2));

    Ain{i} = abs(1./(1+sign(lam1_in(i)-lam1_out(i))*sum(d{i})/2*lam1_in(i)).*1./(1+sign(lam2_in(i)-lam2_out(i))*sum(d{i})/2*lam2_in(i)));
    Aout{i} = abs(1./(1+sign(lam1_out(i)-lam1_in(i))*sum(d{i})/2*lam1_out(i)).*1./(1+sign(lam2_out(i)-lam2_in(i))*sum(d{i})/2*lam2_out(i)));

    alpha = linspace(0,1,layers);
    rho{i} = 1./(Aout{i}-Ain{i}).*(-Ain{i} + sqrt(alpha*Aout{i}.^2+(1-alpha).*Ain{i}.^2));

    rho{i} = sort(rho{i});
    rho{i} = rho{i} - min(rho{i});
    rho{i} = rho{i}./max(rho{i});
    
    D = [0;cumsum(d{i})];
    D = D./max(D);
  
    EVsurf{i}  = interp1(D,fiber{i},rho{i}','linear','extrap');
    EDsurf{i}  = interp1(D,fiber{i},alpha','linear','extrap');
end;


%    EVsurf{FBtri(i,1)}  = interp1(Vx,fiber{FBtri(i,1)}(2:end,:),xi,'pchip','extrap');

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