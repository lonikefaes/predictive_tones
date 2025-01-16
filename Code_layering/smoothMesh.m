function [v2, faces] = smoothMesh(vertices, faces, varargin)
%SMOOTHMESH Smooth mesh by replacing each vertex by the average of its neighbors 
%
%   V2 = smoothMesh(V, F)
%   [V2 F2] = smoothMesh(V, F)
%   Performs smoothing of the values given in V, by using adjacency
%   information given in F. 
%   V is a numeric array representing either vertex coordinate, or value
%   field associated to each vertex. F is an array of faces, given either
%   as a NF-by-3 or NF-by-4 numeric array, or as a cell array. 
%   Artifact adjacencies are added if faces have more than 4 vertices.
%
%   Example
%     [v f] = torusMesh([50 50 50 30 10 30 45]);
%     v = v + randn(size(v));
%     [v2 f] = smoothMesh(v, f, 3);
%     figure; drawMesh(v2, f);
%     l = light; lighting gouraud
%
%   See also
%     meshes3d, meshAdjacencyMatrix, triangulateFaces
%

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2013-04-29,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2013 INRA - Cepia Software Platform.
alpha = 0.4; 
avg = 'laplacian';
avg = 'mean';
% determine number of iterations
nIter = 1;
if ~isempty(varargin)
    nIter = varargin{1};
    if length(varargin)>1
       method = varargin{2};
       
    else
       method = 0;  % 1 for presereving positions at the outer indices
    end
end


if ~method
    % ensure faces correspond to a triangulation
    if iscell(faces) || size(faces, 2) > 3
        faces = triangulateFaces(faces);
    end
    
%    compute adjacency matrix
    adj = meshAdjacencyMatrix(faces);
    
    
    
    
    % Add "self adjacencies"
    nv = size(adj, 1);
    adj = adj + speye(nv);
    
    % weight each vertex by the number of its neighbors
    w = spdiags(full(sum(adj, 2).^(-1)), 0, nv, nv);
    adj = w * adj;
    
    % do averaging to smooth the field
    v2 = vertices;
    for k = 1:nIter
        v2 = adj * v2;
    end
    
else
    %% Old version
    % % Compute vertex adjacencies
    edges = meshEdges(faces);
    FV.vertices = vertices;
    FV.faces = faces;
    Ne   = vertex_neighbours(FV);
    v2   = zeros(size(vertices));
    TR   = triangulation(faces,vertices);
   [~,FBpoints] = freeBoundary(TR);
     FBtri = find(ismember(vertices,FBpoints,'rows')==1);
    % apply several smoothing
    for iter = 1:nIter

        % replace the coords of each vertex by the average coordinate in the
        % neighborhood
        for i = 1:size(vertices, 1),
            N1 = unique(Ne{i});
            N1(N1==i) = [];
%             N12=unique([Ne{Ne{i}}]);
%             N12(N12==i) = [];
            if ~ismember(i,FBtri)
                 switch(avg)
                    case('mean')
                        v2(i,:)  = mean(vertices(N1, :));
                    case('laplacian')
                        v2(i,:) = (1-alpha)*vertices(i,:)+alpha*mean(vertices(N1, :));
                end
            else

                ind =  unique([i,N1(ismember(N1,FBtri))]);
         
                DIFF1 = [diff(vertices(ind,:),[],1)];
                LENG1 = sqrt(sum(DIFF1.^2,2));
                GRAD1 = DIFF1./repmat(LENG1,1,3);
                ANGLE1 =abs(acos(sum(GRAD1(2:end,:).*GRAD1(1:end-1,:),2))).*90./pi;
                %   plot3(vertices(neighInds,1),vertices(neighInds,2),vertices(neighInds,3),'.k');
                %  plot3(vertices(i,1),vertices(i,2),vertices(i,3),'.r');
                if any(abs(ANGLE1)>30) && abs(max(vertices(ind,3)) - min(vertices(ind,3)))>0.2;
                    v2(i, :)  =  vertices(i,:);
                else
                    v2(i, :)  =  [nanmean(vertices(ind,:),1)];
                end
                %  plot3(v2(i,1),v2(i,2),v2(i,3),'.r');
            end
        end
        
        % update for next iteration
        vertices = v2;
        disp(['Iteration no. ',num2str(iter),' finished.'])  
    end
end

if nargout==1
   v2t.vertices = v2;
   v2t.faces = faces;
   clear v2;
   v2 = v2t;
end

return;
     conn = meshconn(faces,size(vertices,1));
    beta  = -1.02*alpha;
    ibeta=1-beta;
    for j=1:iter
        for i=1:nn
            p(idx(i),:)=ialpha*node(idx(i),:)+alpha*mean(node(conn{idx(i)},:)); 
        end
        node=p;
        for i=1:nn
            p(idx(i),:)=ibeta *node(idx(i),:)+beta*mean(node(conn{idx(i)},:)); 
        end
        node=p;
    end

