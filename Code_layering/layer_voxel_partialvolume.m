function [AA LL]= layer_voxel_partialvolume(D,LL,layers,slices)

no_regions   = length(D);
[nx, ny, nz] = size(D{1});
[X Y Z]      = ndgrid(1:nx,1:ny,1:nz);%

AA      = cell(layers-1,1);
[AA{:}] = deal(zeros(nx*ny,nz));
figure;
h0 = tic;
disp('starting');
for k = slices,
    disp(k);
    
    Xg = round(Y(:,:,k));
    Yg = round(X(:,:,k));
    Xg = Xg(:);
    Yg = Yg(:);
    %plot(Xg(:,:)-0.5,Yg(:,:)-0.5,'color',[1  1 1]/2);
    % plot(Xg(:,:)'-0.5,Yg(:,:)'-0.5,'color',[1  1 1]/2);
    for l = 1:layers-1,

        for r = 1:no_regions
            if isstruct(LL{1,1})
                FV     = LL{l,r};
                FV1    = LL{l+1,r};
                try
                    POLYc  = intersectPlaneSurf(FV,[0 0 k+.1],[0 0 1]);
                    
                    POLY1c = intersectPlaneSurf(FV1,[0 0 k+.1],[0 0 1]);
                    
                    
                catch
                    POLYc  = intersectPlaneSurf(FV,[0 0 k-.1],[0 0 1]);
                    
                    POLY1c = intersectPlaneSurf(FV1,[0 0 k-.1],[0 0 1]);
                end
                poly   = POLYc{1}';
                poly1  = POLY1c{1}';
              %  LL{r,l} = POLYc{1}'; 
               % LL{r,l+1} = POLY1c{1}'; 
            else
                poly   = LL{k,l,r};
                poly1  = LL{k,l+1,r};
            end
            
            %clip it
            poly(:,2) = min(poly(:,2),nx); 
            poly(:,2) = max(poly(:,2),1); 
            poly(:,1) = min(poly(:,1),ny); 
            poly(:,1) = max(poly(:,1),1); 
            poly1(:,2) = min(poly1(:,2),nx); 
            poly1(:,2) = max(poly1(:,2),1); 
            poly1(:,1) = min(poly1(:,1),ny); 
            poly1(:,1) = max(poly1(:,1),1); 
            
            poly   = [ones(1,3)*NaN;poly; ones(1,3)*NaN];
            poly1  = [ones(1,3)*NaN;poly1; ones(1,3)*NaN];
            
            [latcells1,loncells1] = polysplit(poly(:,1), poly(:,2));
            [latcells2,loncells2] = polysplit(poly1(:,1),poly1(:,2));
       
            parfor i = 1:length(latcells1);
                lat1 = latcells1{i};
                lon1 = loncells1{i};
                if ~ispolycw(lat1,lon1);
                    [lat1, lon1] = poly2cw(lat1,lon1);
                end;
                lat2 = latcells2{i};
                lon2 = loncells2{i};
                if ~ispolycw(lat2,lon2);
                    [lat2, lon2] = poly2cw(lat2,lon2);
                end;
                latcells{i} = [lat1;lat2(end:-1:1);lat1(1)];
                loncells{i} = [lon1;lon2(end:-1:1);lon1(1)];
            end
            [lat,lon] = polyjoin(latcells,loncells);
            BW        = poly2mask(lat,lon,nx,ny);
       
     
            
            BW(sub2ind([nx,ny],floor(lon),floor(lat))) = 1;
            BW(sub2ind([nx,ny],round(lon),round(lat))) = 1;
            BW(sub2ind([nx,ny],ceil(lon),ceil(lat))) = 1;
            BW = bwmorph(BW,'dilate',1);
            %         imagesc(BW); hold on;
            ind = find(BW==1);
            P = zeros(length(ind),1);
            figure(2), imagesc(BW); colorbar; impixelinfo;
            blocks   = round(linspace(1,length(ind),100));
            h0       = tic;
            
            for p = 1:length(blocks)-1,
                parfor j = blocks(p):blocks(p+1),
                    P(j)= intersect_area(Xg(ind(j)),Yg(ind(j)),lat,lon);
                end
                AA{l}(ind,k) = P;
                h1 = toc(h0);
                disp(['Finished: ',num2str(p+1),'%, ', num2str(h1/60),' minutes elapsed...']);
            end
        end
        
        h1 = toc(h0);
        disp(['Slice = ',num2str(k),', layer = ',num2str(l),', ',num2str(h1/60),' minutes elapsed...']);
        clf;
        tmpA = reshape(AA{l},[nx ny nz]);
        imagesc(tmpA(:,:,k)); colorbar, impixelinfo; hold on
        plot(lat,lon,'g'); drawnow; hold off;
    end
    
end
for l =1:layers-1,
    AA{l} = reshape(AA{l},[nx ny nz]);
end
