function [entire_fiber] = evaluatestreamline3Dsurf(xyz,U_out,VectorF,Roi,b1,b2,opts)



Up   = false;
Down = false;
entire_fiber = [];

entire_fiber = {[],[]};
flag      = 0;
fibers = {[],[]};
for k = 1:2
    % fiber is used to store the coordinates of one fiber
    fiber      = zeros(200,3);
    FU         = zeros(200,3);
    % First fiber coordinate is the current position
    fiber(1,:) = [xyz(1) xyz(2) xyz(3)];
    
    % initialize variables
    check = true; f_length = 1; f_Roi = false; old_gradient = []; f_angle = 0;
    normA_old = Inf;   normB_old = -Inf;
    %  plot3(fiber(1,1),fiber(1,2),fiber(1,3),'.k'); hold on;
    while(check)
        
        gradient  = triinterp8(fiber(f_length,:),VectorF);
        f_U       = triinterp8(fiber(f_length,:),U_out);
        if k ~= 1
            gradient = -gradient;
        end
        gradient  = gradient./(sqrt(sum(gradient.^2))+eps);
        
        % Set a step in the direction of the gradient:
        fiber(f_length+1,:) = fiber(f_length,:) + opts.Step*gradient;
        FU(f_length) = f_U;
        % Determine if this fiber crosses the ROI:
        fiber_r = round(fiber(f_length,:));
        f_Roi   = f_Roi||Roi(fiber_r(1),fiber_r(2),fiber_r(3));
        
        % Check distance to the surface:
        
        % Calculate angle between the current and last on the (fiber) position
        if(f_length>1), f_angle=abs(acos(sum(gradient.*old_gradient)).*180./pi); end;
        
        % Do the Fiber stop checks:
        % Stop if the fiber becomes to long.
        if(f_length>=opts.FiberLengthMax), check=false; %disp('Fiber too long');
        end;
        % Stop if the fiber takes a hard turn.
        if(f_angle>opts.DeviationAngleMax), check=false; %disp('Fiber hard turn');
        end;
        %
        
        % Stop if outside of boundary
        if f_U>(b1*0.99) %||  %|| isnan(f_U)),
            Up = true;
            check=false;
            
        elseif f_U<(b2*1.01)
            Down = true;
            check=false; %disp('Outside of boundary');
        end;
            
            % Check if the fiber will grow outside of the volume
            %  if((sum(fiber(f_length+1,:)>(size(Roi)-1))+sum(fiber(f_length+1,:)<2))>0), check=false;
            %   disp('Outside of ROI');
            %  end;
        
            % Update the fiber length
            if(check), f_length=f_length+1; end;
            
            % Keep the gradient for next step angle change calculation
            old_gradient = gradient;
            
        end,
        
        % Keep the fiber if it is long enough and crossed the Roi
        if((f_length>opts.FiberLengthMin) && f_Roi),
            fibers{k}(1:f_length-1,:) = fiber(1:f_length-1,:);
            fibers{k}(f_length:end,:) = [];
            FUs{k}(1:f_length-1) = FU(1:f_length-1);
            FUs{k}(f_length:end,:) = [];
        else
            fibers{k} = [];
            FUs{k} = [];
        end
        
end

% if (Up==false || Down==false)
%     fibers{1} = []; fibers{2} = [];
% end;


if ~isempty(fibers{1}) && ~isempty(fibers{2})

     for k = 1:2,
         dist   = [sqrt(sum(diff(fibers{k},[],1).^2,2))];
         ds{k} = cumsum([0;dist]);
     end
     
     if ~isempty(ds{2}) && ~isempty(ds{1})
         ED = ds{2}(end)./(ds{1}(end)+ds{2}(end));
     elseif isempty(ds{2})
         ED = ds{1}(end)./(ds{1}(end));
     else
         ED = ds{2}(end)./(ds{2}(end));       
     end
     
     if ED>0.40 && ED<0.6
        if ED~=0.5
            if ED>0.5
                fibers{2} = interp1(ds{2},fibers{2},ds{1},'pchip','extrap');   
            else
                fibers{1} = interp1(ds{1},fibers{1},ds{2},'pchip','extrap');
            end
        end
         
     else
        flag = 1; 
     end

    entire_fiber = fibers;

    if flag
       entire_fiber = {[],[]};
    end

end;
