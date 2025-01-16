function [EV, ED, TC, entire_fiber] = evaluatestreamline3D(xyz,U_out,VectorF,div1,Roi,b1,b2,opts)

EV   = 0;
ED   = 0;
TC   = 0;
Up   = false;
Down = false;
entire_fiber = [];
div  = cell(2,1);
[div{:}] = deal([]);
for k = 1:2
    % fiber is used to store the coordinates of one fiber
    fiber      = zeros(200,3);
    diver      = zeros(200,1);
    % First fiber coordinate is the current position
    fiber(1,:) = [xyz(1) xyz(2) xyz(3)];
    
    % initialize variables
    check = true; f_length = 1; f_Roi = false; old_gradient = []; f_angle = 0;
    normA_old = Inf;   normB_old = -Inf;
    %  plot3(fiber(1,1),fiber(1,2),fiber(1,3),'.k'); hold on;
    while(check)
        
        gradient  = triinterp8(fiber(f_length,:),VectorF);
        divn      = triinterp8(fiber(f_length,:),div1);
        f_U       = triinterp8(fiber(f_length,:),U_out);
        if k ~= 1
            gradient = -gradient;
        end
        gradient  = gradient./(sqrt(sum(gradient.^2))+eps);
        
        % Set a step in the direction of the gradient:
        fiber(f_length+1,:) = fiber(f_length,:) + opts.Step*gradient;
        diver(f_length+1)   = divn;
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
            div{k} = diver(1:f_length-1,:);
        else
            fibers{k} = [];
        end
        
end

% % if length(div{1})>1
% %    div{1} = div{1}(2:end);
% % end
% % if length(div{2})>1
% %    div{2} = div{2}(2:end);
% % end
% divv = [flipud(div{2});div{1}];
% 
% ddivv = divv(round(length(divv)/2)-round(length(divv)/6):round(length(divv)/2)+round(length(divv)/6));

%beta = DCF\divv;
if (Up==false || Down==false)
    fibers{1} = []; fibers{2} = [];
end;

%figure(2)
if ~isempty(fibers{1}) || ~isempty(fibers{2})
%     subplot(121); plot([flipud(div{2}(2:end));div{1}(2:end);]); hold on
%     subplot(122); plot(ddivv); hold on; 
%     
    for k = 1:2,
        dist   = [sqrt(sum(diff(fibers{k},[],1).^2,2))];
        if k == 2
            div{k}  = -div{k};
        end
        S = [cumprod(1 + opts.Step * div{k}, 1)];
        
        cs{k} = cumsum(S);
        ds{k} = cumsum([0;dist]);
    end
    if ~isempty(cs{2}) && ~isempty(cs{1})
        EV = cs{2}(end)./(cs{1}(end)+cs{2}(end));
    elseif isempty(cs{2})
        EV = cs{1}(end)./(cs{1}(end));
    else
        EV = cs{2}(end)./(cs{2}(end));
    end
    if ~isempty(ds{2}) && ~isempty(ds{1})
        ED = ds{2}(end)./(ds{1}(round(length(ds{1})*0.99))+ds{2}(round(length(ds{2})*0.99)));
    elseif isempty(ds{2})
        ED = ds{1}(end)./(ds{1}(end));
    else
        ED = ds{2}(end)./(ds{2}(end));       
    end
    entire_fiber = [fibers{2}(end:-1:2,:);fibers{1}];
    thickness    = [sqrt(sum(diff(entire_fiber,[],1).^2,2))];
    thickness    = cumsum([0;thickness]);
    TC           = max(thickness);
end;
