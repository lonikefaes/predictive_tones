
%%  If you need these codes that implement critical functions with (fast) C code, please visit my website:
%%  http://www.escience.cn/people/gxouyang/Tools.html

%  revise time: May 5 2014, Ouyang,Gaoxiang
%  Email: ouyang@bnu.edu.cn


I = 1;  % RP is the asymmetry matrix

P = rand(1,200);
[RP,DD] = RPplot_FAN(P,3,1,10,0);
subplot(121);imagesc(RP)
[RR1,DET1,ENTR1,L1] = Recu_RQA(RP,I)
P = sin(0.1:0.1:20);
[RP,DD] = RPplot_FAN(P,3,1,10,0);
subplot(122);imagesc(RP)
[RR2,DET2,ENTR2,L2] = Recu_RQA(RP,I)


figure(2)
I = 0;  % RP is the symmetry matrix

P = rand(1,200);
[RP,DD] = RPplot(tp(1:10:end,8)',3,1,.5,0);
subplot(121);imagesc(RP)
[RR1,DET1,ENTR1,L1] = Recu_RQA(RP,I)
P = sin(0.1:0.1:20);
[RP,DD] = RPplot(tp(1:10:end,1)',1,1,.5,0);
subplot(211);imagesc(RP); 
subplot(212);plot([tp(1:10:end-1,1),RP(:,1:10)*20]);
[RR2,DET2,ENTR2,L2] = Recu_RQA(RP,I)
%           norm       : data normalization(0=none, 1=unit interval, 2=zscore, 3=center)
%           Edim       : embedding dimension
%           Tlag       : time lag
%           rescale    : radius rescale option (1=mean distance 2=max distance)
%           radius     : recurrence radius
crp(tp(1:10:end,1));
crqa(tp(1:5:end,1),1,1,.2,80,2)
crqa(tp(1:5:end,1),1,1,.2,80,2,'euc','nonormalize')


figure
N=300; w=40; ws=2;
a=3.4:.6/(N-1):4;
b=.5; for i=2:N, b(i)=a(i)*b(i-1)*(1-b(i-1));end
y=crqa(b,3,2,.1,w,ws);
subplot(2,1,1), plot(a,b,'.','markersize',.1)
title('logistic map'), axis([3.4 4 0 1])
subplot(2,1,2), plot(a(1:ws:N-w),y(1:ws:N-w,1))
ylabel('recurrence rate'), axis([3.4 4 0 1])

aRQA(tp(1:10:end,1), 2, 1, 10, 1, 50, 1, 1);