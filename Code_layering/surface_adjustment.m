function [FVout, fiber] = surface_adjustment(FVin,fibers,Iter)

vertices = FVin.vertices(:,[2 1 3]);
faces = FVin.faces;



alpha= 0.3;
beta = 0.5;

%method = 'mean';
method = 'laplacian';
%method = 'laplacianHC';
%Ne=vertex_neighbours(FVin);
%N=patchnormals(FVin);
for it = 1:Iter
    TR           = triangulation(faces,vertices);
    [~,FBpoints] = freeBoundary(TR);
    FBtri = find(ismember(vertices,FBpoints,'rows')==1);
    v2 = vertices;
    f2 = faces;

    % % Compute vertex adjacencies
    edges = meshEdges(faces);
    
    % replace the coords of each vertex by the average coordinate in the
    
    for i = 1:size(vertices, 1),
      % Nce=unique([Ne{Ne{i}}]);  
        edgeInds  = sum(edges == i, 2) > 0;
        neighInds = unique(edges(edgeInds, :));
        
        if isempty(fibers{i}{1}) || isempty(neighInds)
            v2(i,:) = NaN;
        else
%             if logical(sum(ismember(FBtri,i)))
%  
%                 figure(1),
%                 FVx.vertices = vertices;
%                 FVx.faces = faces;
%                 patch(FVin,'FaceColor',[0 0 1]); view(3); hold on;
%                 plot3(vertices(FBtri,2),vertices(FBtri,1),vertices(FBtri,3),'.g');
%                 plot3(vertices(neighInds,2),vertices(neighInds,1),vertices(neighInds,3),'.r');
%                 
%                 hold off;
%             end
            
            v2(i,:)  = vertices(i,:);
            dist1 = ones(length(neighInds),1)*NaN;
            dist2 = ones(length(neighInds),1)*NaN;
            for k = 1:length(neighInds)
                
                if ~isempty(fibers{neighInds(k)}{1}) %&& i~=neighInds(k),
                    dist1(k) = sum(sqrt(sum(diff(fibers{neighInds(k)}{1},[],1).^2,2)));
                    dist2(k) = sum(sqrt(sum(diff(fibers{neighInds(k)}{2},[],1).^2,2)));
                    if i==neighInds(k),
                        dd1 = [0;cumsum(sqrt(sum(diff(fibers{i}{1},[],1).^2,2)))];
                        dd2 = [0;cumsum(sqrt(sum(diff(fibers{i}{2},[],1).^2,2)))];
                    end
                    %  r2{i}(:,k) = sqrt(sum((fibers{i} - fibers{neighInds(k)}).^2,2));
                end
            end
            
            switch(method)
                case('mean')
                    md1 = nanmean(dist1);
                    md2 = nanmean(dist2);
                case('laplacian')
                    % or laplacian smoothing:
                    md1 = (1-alpha)*dd1(end)+alpha*nanmean(dist1(i~=neighInds));
                    md2 = (1-alpha)*dd2(end)+alpha*nanmean(dist2(i~=neighInds));
%                case('laplacianHC')
                    % or laplacianHC
%                     dm1  = nanmean(dist1(i~=neighInds));
%                     b1   = dm1-(alpha*dd1(end)+(1-alpha)*dm1);
%                     dm1  = dm1-(beta*b1+(1-beta)*nanmean(b1(i~=neighInds)));
%                     dm2  = nanmean(dist2(i~=neighInds));
%                     b2   = dm2-(alpha*dd2(end)+(1-alpha)*dm2);
%                     dm2  = dm2-(beta*b2+(1-beta)*nanmean(b2(i~=neighInds)));
            end
            
            %
            fibers{i}{1}  = interp1(dd1,fibers{i}{1},linspace(0,md1,50),'pchip','extrap');
            fibers{i}{2}  = interp1(dd1,fibers{i}{2},linspace(0,md2,50),'pchip','extrap');
            fiber{i} = [fibers{i}{2}(end:-1:2,:);fibers{i}{1}];
            
        end
    end
    
    indVertices = find(~isnan(v2(:,1)));
    fiber    = {fiber{indVertices}};
    fibers   = {fibers{indVertices}};
    vertices = v2(indVertices,:);
    % create index array for face indices relabeling
    refInds  = zeros(size(indVertices));
    for i = 1:length(indVertices)
        refInds(indVertices(i)) = i;
    end
    
    % select the faces with all vertices within the box
    if isnumeric(f2)
        % Faces given as numeric array
        indFaces = sum(~ismember(f2, indVertices), 2) == 0;
        faces = refInds(f2(indFaces, :));
    end
    
   disp(['Iteration no. ',num2str(it),' finished.'])  
end

FVout.vertices = vertices;
FVout.faces    = faces;
FVout          = MakeContourClockwise3D(FVout);

for i = 1:length(fiber),
   fiber{i} = fiber{i}(:,[2 1 3]); 
end

% FVin  = FV;
% FVout = FV;
% %z 
% FVx  = FV;
% for k = 1:size(fiber{i},1)
%     for i = 1:length(indVertices),
%         ver_in(i,:)  = fiber{i}(k,:);
% %        d(i,:) = sqrt(sum(diff(fiber,[],1).^2,2));
%     end
%     FVx.vertices = ver_in;
%     [Cmean(:,k),Cgaussian(:,k),Dir1,Dir2,Lambda1(:,k),Lambda2(:,k)]=patchcurvature(FVx,true);
%     
% end
% % 
% 
%     ver_in(i,:)  = fiber{i}(1,:);
%     ver_out(i,:) = fiber{i}(end,:);
  %  for i = 1:length(indVertices),
   %     d(i,:) = sqrt(sum(diff(fiber{i},[],1).^2,2));
  %  end
% FVin.vertices = ver_in;
% FVout.vertices = ver_out;
% [Cmean_in,Cgaussian_in,Dir1_in,Dir2_in,Lambda1_in,Lambda2_in]=patchcurvature(FVin,true);
% [Cmean_out,Cgaussian_out,Dir1_out,Dir2_out,Lambda1_out,Lambda2_out]=patchcurvature(FVout,true);

% 
% Lambda1 = abs(Lambda1);
% Lambda2 = abs(Lambda2);
% Ain  = 1./(1+sign(Lambda1(:,1:end-1)-Lambda1(:,2:end)).*d/2.*Lambda1(:,1:end-1)).*(1./(1+sign(Lambda2(:,1:end-1)-Lambda2(:,2:end)).*d/2.*Lambda2(:,1:end-1)));
% Aout = 1./(1+sign(Lambda1(:,2:end)-Lambda1(:,1:end-1)).*d/2.*Lambda1(:,2:end)).*1./(1+sign(Lambda2(:,2:end)-Lambda2(:,1:end-1)).*d/2.*Lambda2(:,2:end));
% 
% 
% Lambda1x = sign(Lambda1).*Lambda1;
% Lambda1x(:,L)
