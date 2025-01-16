function [LS] = get_layer_surfaces(V,layers,M)

%name       = 'name';
no_regions = length(M);
sp         = linspace(0.01,0.99,layers); 

for r = 1:no_regions,
    for i = 1:layers,
        disp(i)
        if no_regions>1
          tmp = V.*M{r} + ~M{r};   
        else
          tmp = V;
        end
        FVx                       = isosurface(tmp,sp(i));
        [FVx.vertices, FVx.faces] = smoothMesh(FVx.vertices, FVx.faces,4,1);
        LS{i,r}                   = FVx;
    end;
end

return
% 
% if nargin<4 
%    display = false; 
%    export = false;
% elseif nargin<5
%    export = false; 
% else
%    export = true;
% end
% 
% 
%    cm         = jet(layers);
% 
% if display, figure; 
%     Data = spm_read_vols(infofile);
%    % Data = Data(:,end:-1:1,:);
%    % Data = Data./max(Data(:));
%    %Data(:,:,68:end)  = [];
%    %Data(:,250:end,:) = [];
%    [nx, ny, nz] = size(Data);
%   % Data(Data==0) = NaN;
%     hz= slice(Data,[],[],[65]);
% %     p1 = patch(isosurface(Data, 0.1),'FaceColor','red',...
% % 	'EdgeColor','none');
% %     p2 = patch(isocaps(Data, 0.1),'FaceColor','interp',...
% % 	'EdgeColor','none');
% %     isonormals(Data,p1)
%     hold on;
%     set(hz,'EdgeColor','none','Facecolor','interp');
%     axis([1 nx 1 ny 1 nz]);
%     %lightangle(0,90)
%     colormap(gray);
%     daspect([1,1,0.4])
%     %axis tight
%     view(113,44);
% %     view( [ -2.4168    4.2717         0   -0.4640;
% %    -4.0640   -2.2993    0.7701    2.6982;
% %    -1.3159   -0.7445   -2.3784   39.1070;
% %          0         0         0    1.0000]);
% %     camlight;
%     camzoom(4.5)
%     camproj perspective;
%     set(hz,'SpecularColorReflectance',0,'SpecularExponent',50)
%  %   caxis([1,2])
% end;     
% %         if display
% %             patch(FVx,'FaceColor',cm(i,:)); view(3);  hold on; %camlight;
% %         end;
%         
%            if export
%             %   g = gifti(LS{i,r});
%             %   gg.cdata = single(repmat(cm(i,:),size(g.vertices,1),1)*255);
%                
%              %  save(g,[name,'_surface_R',num2str(r),'_L',num2str(i),'.vtk'],'Base64Binary');
%                mat           = infofile.mat;
%                MS.vertices   = (mat*[(LS{i,r}.vertices(:,[1 2 3])/(interp_factor/2)+repmat(offset([2 1 3]),size(LS{i,r}.vertices,1),1)),ones(size(LS{i,r}.vertices,1),1)]')';
%              %  MS2.vertices   = ([(LS{i,r}.vertices(:,[2 1 3])/(interp_factor/2)+repmat(offset,size(g.vertices,1),1)),ones(size(g.vertices,1),1)]')';
%                MS.vertices   = MS.vertices(:,1:3); 
%                MS2.vertices   = (LS{i,r}.vertices/(interp_factor/2) + repmat(offset([2 1 3]),size(LS{i,r}.vertices,1),1));
%                %MS2.vertices(:,2)   = ny - MS2.vertices(:,2)+1;
%                MS.faces      = LS{i,r}.faces;
%                MS2.faces     = LS{i,r}.faces;
%                patch(MS2,'FaceColor',cm(i,:),'EdgeColor','k','EdgeAlpha',0); drawnow; %view(3);
%                mvtk_write(MS,[name,'_surface_R',num2str(r),'_L',num2str(i),'.vtp'],'xml-binary');
%               % MS.vector    = repmat(cm(i,:)*255,size(size(LS{i,r}.vertices,1),1),1);
%              %  mvtk_write(MS,[name,'_surface_R',num2str(r),'_L',num2str(i),'.vtk']);
%            end
%     end;
%     
