function show_surface(LS,file,interp_factor,range,slc)

[layers, no_regions] = size(LS);
offset               = [range.r2.mx(1),range.r2.my(1),range.r2.mz(1)];
D                    = load_data(file);
for i = 1:layers
    for r = 1:no_regions
        LS{i,r}.vertices   = (LS{i,r}.vertices/(interp_factor/2) + repmat(offset([2 1 3]),size(LS{i,r}.vertices,1),1));
    end
end
if nargin<5
    slc = mean(LS{1,1}.vertices(:,3));
end


cm  = jet(layers);

figure;
[nx, ny, nz] = size(D.data);

hz = slice(D.data,[],[],[slc]);
hold on;
set(hz,'EdgeColor','none','Facecolor','interp');
axis([1 nx 1 ny 1 nz]);
colormap(gray);
daspect([1,1,0.4]);
view(113,44);
camzoom(4.5)
camproj perspective;
set(hz,'SpecularColorReflectance',0,'SpecularExponent',50)

for i = 1:layers
    for r = 1:no_regions
        patch(LS{i,r},'FaceColor',cm(i,:),'EdgeColor','k','EdgeAlpha',0); drawnow; %view(3);
    end
end

hold off;