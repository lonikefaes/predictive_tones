function [T, XYZtot,XYZup,XYZlow, Width] = laminarity(C,PV,interp_factor,C2,PV2),

layers     = length(PV);
[nx ny nz] = size(PV{1});

T = cell(layers,nz);
XYZtot = cell(layers,nz);
Width  = cell(layers,nz);
for l = 1:layers,
    k = 0;
  
    for sl = 1:nz,
        
        if sum(sum(PV{l}(:,:,sl)))>0,
          
            k   = k+1;
            [xi, yi] = poly2cw(C{sl,l}(:,1),C{sl,l}(:,2));
            [x2i, y2i] = poly2cw(C2{sl,l}(:,1),C2{sl,l}(:,2));
            [x2i1, y2i1] = poly2cw(C2{sl,l+1}(:,1),C2{sl,l+1}(:,2));
            XYZ = equalspacing3d([xi,yi,ones(length(xi),1)]/interp_factor+0.5,20);
            PVi = squeeze(PV{l}(:,:,sl));
            [indXYi] = find(PVi>0);
            [Yi,Xi]  = ind2sub([nx ny],indXYi);
            pos = [];
            DD  = [];
            DD0  = [];
            XYZi = [];
            Dist = [0;cumsum(sqrt(sum(diff(XYZ(:,1:2),[],1).^2,2)))];
            Dist0 = Dist;
            Dist = Dist./max(Dist)*100;
            for i = 1:length(Xi)
                D = sqrt((XYZ(:,1)-Xi(i)).^2+(XYZ(:,2)-Yi(i)).^2);
                [di(i), pos(i)] = min(D);
                DD(i) = Dist(pos(i));
                DD0(i) = Dist0(pos(i));
                XYZi(i,:) = XYZ(pos(i),:);
               % XYZv(i,:)  = [Xi(i),Yi(i),sl];
            end;
            
            [~, ind_sorted] = sort(pos);
            %            figure(2),
            %            plot(XYZ(:,1),XYZ(:,2)); hold on;
            %            plot(Xi,Yi,'.k');
            %            hold off;
            
            figure(1),
            cm = hot(64);
            imagesc(PV{l}(:,:,sl)); colormap(flipud(cm)); hold on; impixelinfo;
            plot(XYZ(:,1),XYZ(:,2));
            for i = 1:length(Xi)
                plot([Xi(i);XYZi(i,1)],[Yi(i);XYZi(i,2)],'k');
            end
       
            

            %   pause;
       
            XYZtot{l,sl}  =  equalspacing3d(XYZ,3000,1);
            UP            =  equalspacing3d([x2i1,y2i1,ones(length(x2i1),1)]/interp_factor+0.5,3000,1);
            LOW           =  equalspacing3d([x2i,y2i,ones(length(x2i),1)]/interp_factor+0.5,3000,1);
            XYZup{l,sl}   =  UP;
            XYZlow{l,sl}  =  LOW;
            
            plot(LOW(:,1),LOW(:,2),'g');
            plot(UP(:,1),UP(:,2),'r');
            plot(XYZtot{l,sl}(1,1),XYZtot{l,sl}(1,2),'oc');
            
            WU = [];
            WL = [];
            for i = 1:size(XYZ,1)
                DL = sqrt((x2i-XYZ(i,1)).^2+(y2i-XYZ(i,2)).^2);
                DU = sqrt((x2i1-XYZ(i,1)).^2+(y2i1-XYZ(i,2)).^2);
                [WL(i)] = min(DL);
                [WU(i)] = min(DU);
            end;
            xi = linspace(1,size(XYZ,1),3000)';
            x  = [1:size(XYZ,1)]';
            Width{l,sl}   = (interp1(x,WL,xi) + interp1(x,WU,xi))/2; 
            T{l,sl} = [Xi(ind_sorted),Yi(ind_sorted),ones(length(Xi),1)*sl,indXYi(ind_sorted),DD(ind_sorted)',PVi(indXYi(ind_sorted)),DD0(ind_sorted)'];
        
            drawnow;
            hold off;
          %  pause;
        end
    end
    %XYZtot{l} = XYZtot{l}/k; 
end