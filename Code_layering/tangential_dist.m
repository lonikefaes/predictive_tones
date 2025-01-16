function Dist = tangential_dist(ALL,M,interp_factor);

[nx,ny,nz] = size(ALL);
no_regions = length(M);
ALL        = logical(ALL);

Mask4               = ones(nx,ny,nz);
Mask4(2:end-1,2:end-1,:) = 0;

MaskN               = ALL.*Mask4;

% figure(1)
% for k = 1:nz
%     subplot(121), imagesc(MaskN(:,:,k)); colorbar; impixelinfo; title('Equivolume')
%     subplot(122), imagesc(M{1}(:,:,k));colorbar; impixelinfo; 
%     title(num2str(k));
%     pause(0.1);
% end
Dist      = cell(no_regions,1);
[Dist{:}] = deal(zeros(nx,ny,nz));
for r = 1:no_regions,
    if no_regions>1
      [EDGE,L]            = bwlabeln(MaskN.*M{r});  % maybe remove .*M{r} for one region
    else
      [EDGE,L]            = bwlabeln(MaskN);     
    end
    D                   = zeros(nx,ny,nz);
    D(EDGE==1)          = 0.0001;
    D(EDGE==2)          = 10000;
    D((ALL.*M{r})==1 & ~MaskN) = -Inf;
    D(D==0)             = NaN;
    if nz>1
    for k=interp_factor+1:nz-interp_factor
        Dist{r}(:,:,k) = laplace2d(D(:,:,k));
    end
    else
        Dist{r}(:,:,1) = laplace2d(D(:,:,1)); 
    end
    Dist{r}(Dist{r}<=0.0001) = 0;
    Dist{r}(Dist{r}>=10000) = 0;
    Dist{r} = Dist{r}./max(Dist{r}(:));
    Dist{r}(~(D==-inf)) = NaN;
end