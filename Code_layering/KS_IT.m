function [SM SP] = KS_IT(DD,m,P,A,H,Q,R,max_iter);       


MM = zeros(size(m,1),size(DD,2));
PP = zeros(size(m,1),size(m,1),size(DD,2));

for it = 1:max_iter
    
    for k=1:size(DD,2)
        [m,P] = kf_predict(m,P,A,Q);
        [m,P] = kf_update(m,P,DD(:,k),H,R);
        
       % [m P F Rx] = vb_glm(DD(:,k),H,0);
      %  R = H*Rx*H';
        MM(:,k) = m;
        PP(:,:,k) = P;
    end
    [SM,SP] = rts_smooth(MM,PP,A,Q);
    m = SM(:,1);
    P = SP(:,:,1);
    disp(it);
end