function FVout = stitch_surfaces(FV1in,FV2in),

vertices1 = FV1in.vertices;
faces1    = FV1in.faces;

vn1       = size(vertices1,1);

vertices2 = FV2in.vertices;
faces2    = FV2in.faces;

vertices  = [vertices1;vertices2];
faces     = [faces1;faces2+vn1];


TR               = triangulation(faces1,vertices1);
[~,FBpoints]     = freeBoundary(TR);

boundary_points  = ismember(vertices1,FBpoints,'rows');
boundary_points  = find(boundary_points==1);
edges = meshEdges(faces1);  
faces3 = [];
faces4 = [];
BP = boundary_points(1);
old_BP = 0;

FVout.vertices = vertices;
FVout.faces = faces;
figure,
patch(FVout,'FaceColor',[1 0 0]); view(3); hold on
f = 0;
for i = 1:length(boundary_points),
    edgeInds  = sum(edges == boundary_points(i), 2) > 0;
    neighInds = unique(edges(edgeInds, :));
    neighInds(neighInds==boundary_points(i)) = [];
    neighbors = neighInds(ismember(neighInds,boundary_points));
    % check for crossing
    %        a = vertices(boundary_points(i),:);
    %        b = vertices(neighbors(1)+vn1,:);
    %        c = vertices(neighbors(1)+vn1,:);
    %        d = vertices(boundary_points(i),:);
    %        theta1 = rad2deg(atan2(norm(cross(a,b)),dot(a,b)));
    %        theta2 = rad2deg(atan2(norm(cross(c,d)),dot(c,d)));
    %if theta1 == theta2
    
%     %        flag = 0;
%     %        if i>1
%     %            if sum(ismember(faces3(:,[1 2]),[boundary_points(i),neighbors(1)],'rows')) || sum(ismember(faces3(:,[1 2]),[neighbors(1), boundary_points(i)],'rows'));
%     %                flag = 1;
%     %                disp(i)
%     %            end
%     %        end
%     %        if flag == 0
%     if i==81;
%        disp(81); 
%     end
%     if old_BP ~= neighbors(1)
%         faces3(i,:) = [BP,neighbors(1),neighbors(1)+vn1];
%         old_BP  = BP;
%         old_neighbors = neighbors;
%         BP      = neighbors(1);
%     else
%         if length(neighbors)>2
%             faces3(i,:) = [BP,neighbors(3),neighbors(3)+vn1];
%             old_BP  = BP;
%             old_neighbors = neighbors;
%             BP      = neighbors(3);
%         else
%             faces3(i,:) = [BP,neighbors(2),neighbors(2)+vn1];
%             old_BP  = BP;
%             old_neighbors = neighbors;
%             BP      = neighbors(2);
%         end
%         
%     end
%     BPx(i) = BP;
    %     faces4(i,:) = [boundary_points(i)+vn1,neighbors(1)+vn1,boundary_points(i)];
    %        else
    if length(neighbors)>2
        cm = {'b','g','m','y'};
        disp(length(neighbors));
        for k = 1:length(neighbors)
        plot3(vertices(neighbors(k),1),vertices(neighbors(k),2),vertices(neighbors(k),3),'.','color',cm{k});
        end
   %    pause; 
    
    else
      f = f+1;  
      faces3(f,:) = [boundary_points(i),neighbors(1),neighbors(1)+vn1];
    end
    %            faces4(i,:) = [boundary_points(i)+vn1,neighbors(2)+vn1,boundary_points(i)];
    %        end
    
    % else
    %  disp('hmmm')
    % end
end

faces = [faces;faces3];
plot3(vertices(neighbors,1),vertices(neighbors,2),vertices(neighbors,3),'.b');

FVout.vertices = vertices;
FVout.faces = faces;
figure,
patch(FVout,'FaceColor',[1 0 0]); view(3); hold on
% 
% for i = 1:length(boundary_points)
%     plot3(vertices(boundary_points(i),1),vertices(boundary_points(i),2),vertices(boundary_points(i),3),'.b');
%   % plot3(vertices(boundary_points(i),1),vertices(boundary_points(i),2),vertices(boundary_points(i),3),'.b');
%   % pause
% end
hold off;
