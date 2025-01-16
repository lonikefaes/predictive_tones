function PV = surface_voxel_partial_volume(FV,D,xi,plot_flag),

[nx ny nz] = size(D);
[xx,yy,zz] = ndgrid(-2:2);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2;

layers  = size(FV,1);

if plot_flag
    figure(1),hold on;
    for i = 1:layers
        patch(FV{i},'FaceColor',[0 0 1]);
    end
end
P = cell(layers,1);
N = cell(layers,1);
[P{:}] = deal(zeros(size(FV{1}.vertices)));
[N{:}] = deal(zeros(size(FV{1}.vertices)));
%     nn2 = fiber{i}(1:end-1,:) - fiber{i}(2:end,:);
%     nn2 = nn2./repmat(sqrt(sum(nn2.^2,2)),1,3);
%     NN2{i} = interp1([1:size(nn2)]',nn2,[0:size(nn2)]','linear','extrap');

for i = 1:layers
    P{i}  = FV{i}.vertices;
   % N{i}  = -patchnormals(FV{i});   
   if i<layers
       N{i} = FV{i}.vertices - FV{i+1}.vertices;
       N{i} = N{i}./repmat(sqrt(sum(N{i}.^2,2)),1,3);
   else
       N{i} = N{i-1};
   end
    if plot_flag
       % if i==1
            h = quiver3(P{i}(:,1), P{i}(:,2), P{i}(:,3),N{i}(:,1), N{i}(:,2), N{i}(:,3)); hold on;
       % end
    end
end;

P       = cell2mat(P);
N       = cell2mat(N);
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
parfor i = 1:length(L) 
    L{i}   = zeros(1,layers-1);
end
PVv      = L;
M        = kron(eye(layers),ones(1,3)/3);
for p = 1:length(blocks)-1,
    parfor k = blocks(p):blocks(p+1),
        [L{k}, PVv{k}] = evaluate_partial_volume(FV,P,N,PP,[x0(k),y0(k),z0(k)],f,e,layers,M,plot_flag);
    end;
    h1 = toc(h0);
    disp(['Finished: ',num2str(p+1),'%, ', num2str(h1/60),' minutes elapsed...']);
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
XX =  zeros(nx*ny*nz,layers-1);
for i = 1:layers-1
  XX(:,i)  = PV{i}(:);
end;
ind   = find(sum(XX,2)>=1);
XX(ind,:) = XX(ind,:)./repmat(sum(XX(ind,:),2),1,layers-1);
XX(XX>1) = 1;
ind   = find(sum(XX,2)<1 & sum(XX(:,[1,layers-1]),2)==0);
XX(ind,:) = XX(ind,:)./repmat(sum(XX(ind,:),2),1,layers-1);

for i = 1:layers-1
  PV{i}(:)  = XX(:,i);
end;


end % function
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function con = condition(vertices,faces)

con = ~isempty(vertices) && ~isempty(faces);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, Volume] = evaluate_partial_volume(FV,P,N,PP,xyz,f,e,layers,M,plot_flag),

L      = zeros(1,layers-1);
Volume = zeros(1,layers-1);

%if ~isnan(d)
    % define cube:
    [v, box]  = cube(xyz(1),xyz(2),xyz(3),0);
    mv        = mean(v,1);

    [sd, d, points,normals]   = signed_distance(P,mv,N,PP,layers,M);

    % compute the edge list
    cube_edges  = [ v(e(:,1), :) v(e(:,2), :) ];

    ind        = find(abs(sd)<=2);
    [~,ind0]   = sort(sd(ind),'descend');
    ind        = ind(ind0);
    if plot_flag
        fv.vertices = v; fv.faces = f;
        figure(2); clf; set(gcf, 'renderer', 'opengl');
        patch(fv,'FaceColor',[1 1 0],'FaceAlpha',0.2); view(3); hold on;
    end

    Volume_parts = zeros(1,layers);
    %PLANE = [];
    col  = jet(layers);
    for i = ind(:)',%1:layers
        
        plane = createPlane(points(i,:),normals(i,:));
        %PLANE = [PLANE;plane]; 
        %axis([mv(1)-3 mv(1)+3 mv(2)-3 mv(2)+3 mv(3)-3 mv(3)+3]);
       % drawPlane3d(plane,col(i,:)); drawnow
 
        % identify which edges cross the mesh
        inds = isBelowPlane(v, plane);
        IntPoints = [];
        if sum(inds)<8 && sum(inds)~=0
            edgeCrossInds = find(sum(inds(e), 2) == 1);
            % compute one intersection point for each edge
            IntPoints = intersectEdgePlane(cube_edges(edgeCrossInds, :), plane);
        else
            IntPoints = [];
        end;
        if ~isempty(IntPoints)
            [K, Volume_Below] = convhulln([v(inds,:);IntPoints]);
            Volume_parts(i)   = Volume_Below;
        end
    end
    if plot_flag
        hold off;
        %     figure(5),
    end
    Volume_diff = zeros(1,layers);
    if sum(Volume_parts)~=0,
        st = find(Volume_parts~=0);
        Volume_diff(st(1):end) = abs(diff([1,Volume_parts(st(1):end)]));
        Volume(:)              = Volume_diff(2:end) ;
%         if sum(Volume)>1 || (sum(Volume)<1 && sum(Volume([1,layers-1]))==0)
%            disp(Volume); 
%         end;
%         
%        % disp(Volume_parts);s
%         if sum(Volume)==1
%         disp(sum(Volume));
%         else
%             disp([sum(Volume),st]);  
%         end
        
        L(find(Volume>0))      = find(Volume>0);
 
    else
        
        sd         = logical([diff(sign(sd))]);
        L(sd)      = find(sd==1);
        Volume(sd) = 1;

    end;
      %  disp(Volume);
      %  disp(L);
%figure(3), plot(mm)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signed_d, d,p,n] = signed_distance(P,mv,N,PP,layers,M),
    PP(:)     = sqrt((P(:,1)-mv(1)).^2 + (P(:,2)-mv(2)).^2 + (P(:,3)-mv(3)).^2);
    [d,I]     = sort(PP,1);
    I         = I(1,:);
    I         = I+[0:size(PP,1):(layers-1)*size(PP,1)];%repmat(,1,1);
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