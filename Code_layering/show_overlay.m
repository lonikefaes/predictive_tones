function show_overlay(D0,LL, interp_factor,cm,mask),

figure; 
clf;
dim = size(D0);

if length(dim)<4
    dim(4) = 1;
end

if nargin<5
    M = [];
else
    ind = find(mask == 1);
    [M(:,1),M(:,2),M(:,3)] = ind2sub([dim(1:3)],ind);
    M(:,1:2) = M(:,1:2) - 0.5;
end

     D  = D0;
    tx  = 0;
    slices =  logical(squeeze(sum(sum(mean(D,4),1),2)));
    [ha, pos] = tight_subplot(sum(slices),dim(4),[.001 .001],[.01 .001],[.001 .001]);
 
    nzi = 1:dim(3);
    for k = nzi(slices),
  
        for t = 1:dim(4)
            tx = tx + 1;
            axes(ha(tx));
            imagesc(D(:,:,k,t),[-4 10]); hold on; %axis image;
            colormap(eval(cm)); impixelinfo;
            if k>1 && k<=dim(3)-1,
                % POLY = intersectPlaneSurf(LS{1},[0 0 k*interp_factor],[0 0 1]);
                plot(LL{k,end}(:,1)/interp_factor,LL{k,end}(:,2)/interp_factor,'w','linewidth',2);
                % POLY = intersectPlaneSurf(LS{end},[0 0 k*interp_factor],[0 0 1]);
                plot(LL{k,1}(:,1)/interp_factor,LL{k,1}(:,2)/interp_factor,'w','linewidth',2);
                if any(M(:,3)==k)
                    ind = (M(:,3)==k);
                    plot(M(ind,2),M(ind,1),'ok','linewidth',1,'MarkerSize',5);
                end
            end
            hold off;
            set(ha(tx),'XTickLabel','','YTickLabel','');
            
            %  end
        end
        if tx == sum(slices)*dim(4);
            h    = colorbar('location','south');
            pos  = get(ha(tx),'Position');
            set(h, 'Position', [pos(1)*0.8 pos(2)+0.05 0.2 pos(4)*0.2])
        end
    end