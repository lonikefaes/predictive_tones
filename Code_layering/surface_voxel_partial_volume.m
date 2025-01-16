function PV = surface_voxel_partial_volume(FV,D,xi,plot_flag),

[nx ny nz] = size(D);
[xx,yy,zz] = ndgrid(-2:2);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1;

layers  = size(FV,1);
P      = [];
N      = [];
for i = 1:layers
    P  = [P;FV{i}.vertices];
    N  = [N;patchnormals(FV{i})];   
end;
PP      = zeros(size(FV{1}.vertices,1),layers);

PV      = cell(layers-1,1);
[PV{:}] = deal(zeros(nx,ny,nz));

Mask = imdilate(D>=xi(1) & D<=xi(end),nhood);

IndM       = find(Mask==1);
[y0,x0,z0] = ind2sub([nx,ny,nz],IndM);

% create voxel as a mesh
[v, f]   = createCube;
e        = meshEdges(f);

blocks   = round(linspace(1,length(IndM),100));
h0       = tic;
L        = cell(length(IndM),1);
PVv      = cell(length(IndM),1);
[L{:}]   = deal(zeros(1,layers));
[PVv{:}] = deal(zeros(1,layers));

for p = 1:length(blocks)-1,
    parfor k = blocks(p):blocks(p+1),
        [L{k}, PVv{k}] = evaluate_partial_volume(FV,P,N,PP,[x0(k),y0(k),z0(k)],f,e,layers,plot_flag);
    end
    
    h1 = toc(h0);
    disp(['Finished: ',num2str(p),'%, ', num2str(h1/60),' minutes elapsed...']);
end

for k = 1:length(IndM)
    ind = find(L{k}~=0);
    if ~isempty(ind)
        for j = 1:length(ind)
            if j<layers
                PV{L{k}(ind(j))}(IndM(k)) = PVv{k}(ind(j));
            end
        end
    end
end

end % function
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function con = condition(vertices,faces)

con = ~isempty(vertices) && ~isempty(faces);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, Volume] = evaluate_partial_volume(FV,P,N,PP,xyz,f,e,layers,plot_flag),

L      = zeros(1,layers);
Volume = zeros(1,layers);

%if ~isnan(d)
    % define cube:
    [v, box]  = cube(xyz(1),xyz(2),xyz(3),0);
    mv        = mean(v,1);

    [sd, d] = signed_distance(P,mv,N,PP,layers);

    % compute the edge list
    cube_edges  = [ v(e(:,1), :) v(e(:,2), :) ];
    
    % choose  closest surfaces;
    [~, Sind]  = sort(d);
    Sind       = sort(Sind(1:ceil(layers/2)),'descend');
    FVa        = {};
    ind        = [];
    if plot_flag
        fv.vertices = v; fv.faces = f;
        figure(1); clf; set(gcf, 'renderer', 'opengl');
        patch(fv,'FaceColor',[1 1 0],'FaceAlpha',0.2); view(3); hold on;
    end
    for i = 1:length(Sind)
        vertices   = FV{Sind(i)}.vertices;
        faces      = FV{Sind(i)}.faces;  % faces are same for all;
        
        % 1:length(Vind)
        %     Dx  = double(Mask);
        %     Dx(y0(k),x0(k),z0(k)) = 2;
        
        
        %  slice(Dx,[],[],z0(k)); hold on;
        
        [FVa{i}.vertices, FVa{i}.faces] = clipMeshVerticesSphere(vertices, faces, mean(v,1),2);
        if ~isempty(FVa{i}.vertices),
            ind = [ind,i];
        end
    end
    Volume_parts = [];
    for i = ind
        if ~isempty(FVa{i}.vertices) && plot_flag
            color = ones(1,3);
            if i<=3
                color(i) = 0;
            else
                color(6-i) = 0;
            end
            patch(FVa{i},'FaceColor',color);
        end
        
        if condition(FVa{i}.vertices,FVa{i}.faces)
            if size(FVa{i}.faces,1)==3 && size(FVa{i}.faces,2)==1
                FVa{i}.faces = FVa{i}.faces';
            end
            
            [IntPoints,TriPoints,BlwPoints] = get_points(FVa{i}.faces,FVa{i}.vertices,f,v,e,cube_edges,box,plot_flag);
            VolumeBelow                     = get_volume(IntPoints,TriPoints,BlwPoints,i,plot_flag);
            Volume_parts                    = [Volume_parts;VolumeBelow];
        else
            Volume_parts = [Volume_parts;0];
        end;
    end
    if plot_flag
        hold off;
        %     figure(5),
    end
    
    st = find(Volume_parts ~= 0,1);
    if ~isempty(st),
        L(Sind(ind(st:end)))      = Sind(ind(st:end));
        Volume(Sind(ind(st:end))) = abs(diff([1;Volume_parts(st:end)]));
        L(L==layers)              = 0;
        Volume(L==layers)         = 0;
     
        %             for i = 1:length(Sind)
        %                 PV{L(i)}(y0(k),x0(k),z0(k)) = Volume(i);
        %                 %    imagesc(PV{Sind(i)}(:,:,z0(k))); colorbar; impixelinfo; drawnow;
        %             end
    else
        
        sd         = logical([diff(sign(sd));0]);
        L(sd)      = find(sd==1);
        Volume(sd) = 1;

    end
%figure(3), plot(mm)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signed_d, d] = signed_distance(P,mv,N,PP,layers),
    PP(:)     = sqrt((P(:,1)-mv(1)).^2 + (P(:,2)-mv(2)).^2 + (P(:,3)-mv(3)).^2);
    [d,I]     = min(PP,[],1);
    I         = I+[0:size(PP,1):(layers-1)*size(PP,1)];
    p         = P(I,:);
    n         = N(I,:);
    % Uses Hessian form, ie : N.p = d
    % I this case, d can be found as : -N.p0, when N is normalized
    signed_d = -sum(bsxfun(@times, n, bsxfun(@minus,p,mv)), 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IntPoints,TriPoints,BlwPoints] = get_points(faces,vertices,f,v,e,cube_edges,box,plot_flag),

IntPoints   = [];
TriPoints   = [];
BlwPoints   = [];
Planes      = [];
for i = 1:size(faces,1)
    
    Tri   = vertices(faces(i,:),:);
   
    plane = createPlane(Tri(1,:),Tri(2,:),Tri(3,:));
   
    % identify which edges cross the mesh
    inds = isBelowPlane(v, plane);
    
    if sum(inds)<8
        edgeCrossInds = find(sum(inds(e), 2) == 1);
        % compute one intersection point for each edge
        IntersectionPoints = intersectEdgePlane(cube_edges(edgeCrossInds, :), plane);
    else
        IntersectionPoints = [];
    end;
    
    if ~isempty(IntersectionPoints),
        
        out = isPointInTriangle(IntersectionPoints,vertices,faces(i,:));
        
        if any(out==1)
            Planes     = [Planes;plane];
            IntPoints  = [IntPoints; IntersectionPoints(out==1,:)];
            TriPoints  = [TriPoints; Tri];
        end
    end
end


if ~isempty(IntPoints)
    

    mTriPoints  = mean(TriPoints,1);
    mPlane      = mean(Planes,1);
    mPlane(1:3) = mTriPoints;

    inds       = isBelowPlane(v,mPlane);
    edgeCrossInds = find(sum(inds(e), 2) == 1);
    IntPoints = [IntPoints; intersectEdgePlane(cube_edges(edgeCrossInds, :),mPlane)];
    
    BlwPoints = v(logical(inds),:);

    TriPoints = clipPoints3d(TriPoints,box);
    if plot_flag
        drawPlane3d(mPlane,[1 0 1]);
        plot3(IntPoints(:,1),IntPoints(:,2),IntPoints(:,3),'.m');
        plot3(TriPoints(:,1),TriPoints(:,2),TriPoints(:,3),'.k');
        plot3(BlwPoints(:,1),BlwPoints(:,2),BlwPoints(:,3),'.r');
        
    end
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = get_volume(IntPoint,TriPoints,BlwPoints,c,plot_flag),

Points = [unique(IntPoint,'rows');unique(TriPoints,'rows');unique(BlwPoints,'rows')];
if size(Points,1)>3
    [K, V] = convhulln(Points);
    FVb.vertices = Points;
    FVb.faces    = K;
    if plot_flag
    color = zeros(1,3);
    if c<3
        color(c) = 1;
    else
        color(6-c) = 1;
    end
    patch(FVb,'FaceColor',color,'FaceAlpha',0.7);
    end
else
    V = 0;
end

end % function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [points, lims] = cube(x,y,z,offset),


points = [x y z; ...
          x+1 y z; ...
          x y+1 z; ...
          x+1 y+1 z; ...
          x y z+1; ...
          x+1 y z+1; ...
          x y+1 z+1; ...
          x+1 y+1 z+1];

lims  = [min(points(:,1))-offset,max(points(:,1))+offset,...
         min(points(:,2))-offset,max(points(:,2))+offset,...
         min(points(:,3))-offset,max(points(:,3))+offset];

end