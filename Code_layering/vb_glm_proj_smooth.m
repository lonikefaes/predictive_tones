function [ws Cs F w C] = vb_glm_proj_smooth(y,X,W,ard)


T = size(y,2);

for t = 1:T
    
    if t==1
        [w(:,t) C(:,:,t) F(t) R(:,:,t) a0 b0 c0 d0] = vb_glm_proj(y(:,t),X,W,ard);
    else
        [w(:,t) C(:,:,t) F(t) R(:,:,t) a0 b0 c0 d0] = vb_glm_proj(y(:,t),X,W,ard,a0,b0,c0,d0);
    end
end



[ws,Cs] = rts_smooth(w,C,eye(size(C,1)),mean(R,3)/5);

F = sum(F);
