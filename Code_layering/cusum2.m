function W = cusum2(Y,meanY,tdev,mshift,show)

for j = 1:length(Y)
    
    
    for r = 1:2,
        if r==1
            y = Y(j:-1:1);
        else
            y = Y(j:end);
        end
        tmean = meanY(j);
        
        uppersum = zeros(size(y));
        lowersum = zeros(size(y));
        
        uppersum(1) = 0;
        lowersum(1) = 0;
   %     mshift = 2;
        climit = 1;
        
        for i=2:length(y)
            uppersum(i) = max(0, uppersum(i-1) + y(i) - tmean - mshift*tdev/2);
            lowersum(i) = min(0, lowersum(i-1) + y(i) - tmean + mshift*tdev/2);
        end
        
        
        iupper = find(uppersum >  climit*tdev, 1, 'first');
        ilower = find(lowersum < -climit*tdev, 1, 'first');
 %       min(iupper,ilower)
        if isempty(iupper), iupper = length(y); end
        if isempty(ilower), ilower = length(y); end
        if r~=1
       %         figure(f1),
       % subplot(2,1,r); plot([uppersum,lowersum]); 
               Window(r) = j + min(iupper,ilower)-1;
        else
            
               Window(r) = j - min(iupper,ilower)+1; 
            %       figure(f1),
       % subplot(2,1,r); plot([uppersum,lowersum]);
        end
     %   drawnow;
    end
    Window = sort(Window);
    if j>1
    if Window(1)<=Window_old(1)
       Window(1) = Window_old(1);
    end
    end
    W(j,1) = Window(2)- Window(1);
    Window_old = Window;
    
    if show
        figure(2),
        plot([Y(:),meanY(:)]); hold on;
        plot(j,Y(j),'or')
        w = zeros(size(Y));
        w(Window(1):Window(2)) = 10;
        plot(w,'g'); hold off
        drawnow;
    end
  
end

W =  W./length(W);



