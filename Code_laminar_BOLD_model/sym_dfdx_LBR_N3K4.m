function dfdx_s = sym_dfdx_LBR_N3K4(xn1_1,xn2_1,xn3_1,xn1_2,xn2_2,xn3_2,xn1_3,xn2_3,xn3_3,xn1_4,xn2_4,xn3_4,xk1_1,xk2_1,xk3_1,xk4_1,xk1_2,xk2_2,xk3_2,xk4_2,xk1_3,xk2_3,xk3_3,xk4_3,xk1_4,xk2_4,xk3_4,xk4_4,a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3,mu1,mu2,mu3,cu1,cu2,cu3,lam1,lam2,lam3,c11,c12,c13,c21,c22,c23,c31,c32,c33,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k1_3,n2k2_3,n2k3_3,n2k4_3,al_v1,al_v2,al_v3,al_v4,al_d1,al_d2,al_d3,al_d4,tt_v1,tt_v2,tt_v3,tt_v4,tt_d1,tt_d2,tt_d3,tt_d4,nr1,nr2,nr3,nr4,ve_v1,ve_v2,ve_v3,ve_v4,ve_d1,ve_d2,ve_d3,ve_d4,f0v1,f0v2,f0v3,f0v4,f0d1,f0d2,f0d3,f0d4)
%SYM_DFDX_LBR_N3K4
%    DFDX_S = SYM_DFDX_LBR_N3K4(XN1_1,XN2_1,XN3_1,XN1_2,XN2_2,XN3_2,XN1_3,XN2_3,XN3_3,XN1_4,XN2_4,XN3_4,XK1_1,XK2_1,XK3_1,XK4_1,XK1_2,XK2_2,XK3_2,XK4_2,XK1_3,XK2_3,XK3_3,XK4_3,XK1_4,XK2_4,XK3_4,XK4_4,A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,MU1,MU2,MU3,CU1,CU2,CU3,LAM1,LAM2,LAM3,C11,C12,C13,C21,C22,C23,C31,C32,C33,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K1_3,N2K2_3,N2K3_3,N2K4_3,AL_V1,AL_V2,AL_V3,AL_V4,AL_D1,AL_D2,AL_D3,AL_D4,TT_V1,TT_V2,TT_V3,TT_V4,TT_D1,TT_D2,TT_D3,TT_D4,NR1,NR2,NR3,NR4,VE_V1,VE_V2,VE_V3,VE_V4,VE_D1,VE_D2,VE_D3,VE_D4,F0V1,F0V2,F0V3,F0V4,F0D1,F0D2,F0D3,F0D4)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    10-Mar-2020 11:01:33

t2 = 1.0./xn1_4;
t3 = 1.0./xn2_4;
t4 = 1.0./xn3_4;
t5 = 1.0./tt_v1;
t6 = 1.0./xk1_1;
t7 = tt_v1+ve_v1;
t8 = 1.0./t7;
t9 = xn1_4-1.0;
t10 = xn2_4-1.0;
t11 = xn3_4-1.0;
t12 = n2k1_1.*t9;
t13 = n2k1_2.*t10;
t14 = n2k1_3.*t11;
t15 = 1.0./al_v1;
t16 = 1.0./tt_v2;
t17 = 1.0./xk2_1;
t18 = tt_v2+ve_v2;
t19 = 1.0./t18;
t20 = n2k2_1.*t9;
t21 = n2k2_2.*t10;
t22 = n2k2_3.*t11;
t23 = 1.0./al_v2;
t24 = 1.0./tt_v3;
t25 = 1.0./xk3_1;
t26 = tt_v3+ve_v3;
t27 = 1.0./t26;
t28 = n2k3_1.*t9;
t29 = n2k3_2.*t10;
t30 = n2k3_3.*t11;
t31 = 1.0./al_v3;
t32 = 1.0./tt_v4;
t33 = 1.0./xk4_1;
t34 = tt_v4+ve_v4;
t35 = 1.0./t34;
t36 = n2k4_1.*t9;
t37 = n2k4_2.*t10;
t38 = n2k4_3.*t11;
t39 = 1.0./al_v4;
t40 = 1.0./xk1_2;
t41 = 1.0./nr1;
t42 = 1.0./xk1_1.^2;
t43 = t12+t13+t14+1.0;
t44 = t43.*ve_v1;
t45 = xk1_1.^t15;
t46 = t45.*tt_v1;
t47 = t44+t46;
t48 = t15-1.0;
t49 = xk1_1.^t48;
t50 = 1.0./xk2_2;
t51 = 1.0./nr2;
t52 = 1.0./xk2_1.^2;
t53 = t20+t21+t22+1.0;
t54 = t53.*ve_v2;
t55 = xk2_1.^t23;
t56 = t55.*tt_v2;
t57 = t54+t56;
t58 = t23-1.0;
t59 = xk2_1.^t58;
t60 = 1.0./xk3_2;
t61 = 1.0./nr3;
t62 = 1.0./xk3_1.^2;
t63 = t28+t29+t30+1.0;
t64 = t63.*ve_v3;
t65 = xk3_1.^t31;
t66 = t65.*tt_v3;
t67 = t64+t66;
t68 = t31-1.0;
t69 = xk3_1.^t68;
t70 = 1.0./xk4_2;
t71 = 1.0./nr4;
t72 = 1.0./xk4_1.^2;
t73 = t36+t37+t38+1.0;
t74 = t73.*ve_v4;
t75 = xk4_1.^t39;
t76 = t75.*tt_v4;
t77 = t74+t76;
t78 = t39-1.0;
t79 = xk4_1.^t78;
t80 = 1.0./f0d1;
t81 = 1.0./f0d2;
t82 = 1.0./f0d3;
t83 = f0v1.*n2k1_1.*t8.*t80.*ve_v1;
t84 = tt_d2+ve_d2;
t85 = 1.0./t84;
t86 = f0v2.*n2k2_1.*t19.*t81.*ve_v2;
t87 = tt_d3+ve_d3;
t88 = 1.0./t87;
t89 = f0v3.*n2k3_1.*t27.*t82.*ve_v3;
t90 = tt_d4+ve_d4;
t91 = 1.0./t90;
t92 = f0d4.*n2k4_1.*t35.*t82.*t91.*ve_d4.*ve_v4;
t93 = t89+t92;
t94 = f0d3.*t81.*t88.*t93.*ve_d3;
t95 = t86+t94;
t96 = f0d2.*t80.*t85.*t95.*ve_d2;
t97 = 1.0./tt_d1;
t98 = 1.0./xk1_3;
t99 = tt_d1+ve_d1;
t100 = 1.0./t99;
t101 = f0v1.*n2k1_2.*t8.*t80.*ve_v1;
t102 = f0v2.*n2k2_2.*t19.*t81.*ve_v2;
t103 = f0v3.*n2k3_2.*t27.*t82.*ve_v3;
t104 = f0d4.*n2k4_2.*t35.*t82.*t91.*ve_d4.*ve_v4;
t105 = t103+t104;
t106 = f0d3.*t81.*t88.*t105.*ve_d3;
t107 = t102+t106;
t108 = f0d2.*t80.*t85.*t107.*ve_d2;
t109 = f0v1.*n2k1_3.*t8.*t80.*ve_v1;
t110 = f0v2.*n2k2_3.*t19.*t81.*ve_v2;
t111 = f0v3.*n2k3_3.*t27.*t82.*ve_v3;
t112 = f0d4.*n2k4_3.*t35.*t82.*t91.*ve_d4.*ve_v4;
t113 = t111+t112;
t114 = f0d3.*t81.*t88.*t113.*ve_d3;
t115 = t110+t114;
t116 = f0d2.*t80.*t85.*t115.*ve_d2;
t117 = f0v1.*t8.*t47.*t80;
t118 = f0v2.*t19.*t57.*t81;
t119 = 1.0./al_d3;
t120 = xk3_3.^t119;
t121 = t120.*tt_d3;
t122 = f0v3.*t27.*t67.*t82;
t123 = 1.0./al_d4;
t124 = xk4_3.^t123;
t125 = t124.*tt_d4;
t126 = t35.*t77.*ve_d4;
t127 = t125+t126;
t128 = f0d4.*t82.*t91.*t127;
t129 = t122+t128;
t130 = t129.*ve_d3;
t131 = t121+t130;
t132 = f0d3.*t81.*t88.*t131;
t133 = t118+t132;
t134 = t133.*ve_d2;
t135 = 1.0./al_d2;
t136 = xk2_3.^t135;
t137 = t136.*tt_d2;
t138 = t134+t137;
t139 = f0d2.*t80.*t85.*t138;
t140 = 1.0./al_d1;
t141 = t135-1.0;
t142 = xk2_3.^t141;
t143 = t119-1.0;
t144 = xk3_3.^t143;
t145 = t123-1.0;
t146 = xk4_3.^t145;
t147 = 1.0./tt_d2;
t148 = 1.0./xk2_3;
t149 = 1.0./tt_d3;
t150 = 1.0./xk3_3;
t151 = 1.0./tt_d4;
t152 = 1.0./xk4_3;
t153 = 1.0./f0d4;
t154 = t83+t96;
t155 = 1.0./xk1_4;
t156 = t101+t108;
t157 = t109+t116;
t158 = 1.0./xk1_3.^2;
t159 = xk1_3.^t140;
t160 = t159.*tt_d1;
t161 = t117+t139;
t162 = t161.*ve_d1;
t163 = t160+t162;
t164 = t140-1.0;
t165 = xk1_3.^t164;
t166 = 1.0./xk2_3.^2;
t167 = 1.0./xk2_4;
t168 = 1.0./xk3_3.^2;
t169 = 1.0./xk3_4;
t170 = 1.0./xk4_3.^2;
t171 = 1.0./xk4_4;
t172 = t35.*t72.*t77.*xk4_2;
dfdx_s = reshape([a1_1,a2_1,a3_1,lam1,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,a1_2,a2_2,a3_2,0.0,lam2,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,a1_3,a2_3,a3_3,0.0,0.0,lam3,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu1,0.0,0.0,-lam1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu2,0.0,0.0,-lam2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu3,0.0,0.0,-lam3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c11,0.0,0.0,c21.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c12,0.0,0.0,c22.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c13,0.0,0.0,c23.*t4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn1_4.^2.*(c31.*t9-c21.*xn1_3)-c31.*t2,0.0,0.0,t5.*t6.*(n2k1_1-n2k1_1.*t8.*ve_v1),t16.*t17.*(n2k2_1-n2k2_1.*t19.*ve_v2),t24.*t25.*(n2k3_1-n2k3_1.*t27.*ve_v3),t32.*t33.*(n2k4_1-n2k4_1.*t35.*ve_v4),t5.*t40.*(n2k1_1.*t41-n2k1_1.*t6.*t8.*ve_v1.*xk1_2),t16.*t50.*(n2k2_1.*t51-n2k2_1.*t17.*t19.*ve_v2.*xk2_2),t24.*t60.*(n2k3_1.*t61-n2k3_1.*t25.*t27.*ve_v3.*xk3_2),t32.*t70.*(n2k4_1.*t71-n2k4_1.*t33.*t35.*ve_v4.*xk4_2),t97.*t98.*(t83+t96-t100.*t154.*ve_d1),t147.*t148.*(t86+t94-t85.*t95.*ve_d2),t149.*t150.*(t89+t92-t88.*t93.*ve_d3),t151.*t152.*(f0v4.*n2k4_1.*t35.*t153.*ve_v4-n2k4_1.*t35.*t91.*ve_d4.*ve_v4),t97.*t155.*(-t98.*t100.*t154.*ve_d1.*xk1_4+f0v1.*n2k1_1.*t6.*t8.*t80.*ve_v1.*xk1_2+f0d2.*t80.*t85.*t95.*t148.*ve_d2.*xk2_4),t147.*t167.*(-t85.*t95.*t148.*ve_d2.*xk2_4+f0v2.*n2k2_1.*t17.*t19.*t81.*ve_v2.*xk2_2+f0d3.*t81.*t88.*t93.*t150.*ve_d3.*xk3_4),t149.*t169.*(-t88.*t93.*t150.*ve_d3.*xk3_4+f0v3.*n2k3_1.*t25.*t27.*t82.*ve_v3.*xk3_2+f0d4.*n2k4_1.*t35.*t82.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),t151.*t171.*(n2k4_1.*t33.*t35.*ve_v4.*xk4_2-n2k4_1.*t35.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn2_4.^2.*(c32.*t10-c22.*xn2_3)-c32.*t3,0.0,t5.*t6.*(n2k1_2-n2k1_2.*t8.*ve_v1),t16.*t17.*(n2k2_2-n2k2_2.*t19.*ve_v2),t24.*t25.*(n2k3_2-n2k3_2.*t27.*ve_v3),t32.*t33.*(n2k4_2-n2k4_2.*t35.*ve_v4),t5.*t40.*(n2k1_2.*t41-n2k1_2.*t6.*t8.*ve_v1.*xk1_2),t16.*t50.*(n2k2_2.*t51-n2k2_2.*t17.*t19.*ve_v2.*xk2_2),t24.*t60.*(n2k3_2.*t61-n2k3_2.*t25.*t27.*ve_v3.*xk3_2),t32.*t70.*(n2k4_2.*t71-n2k4_2.*t33.*t35.*ve_v4.*xk4_2),t97.*t98.*(t101+t108-t100.*t156.*ve_d1),t147.*t148.*(t102+t106-t85.*t107.*ve_d2),t149.*t150.*(t103+t104-t88.*t105.*ve_d3),t151.*t152.*(f0v4.*n2k4_2.*t35.*t153.*ve_v4-n2k4_2.*t35.*t91.*ve_d4.*ve_v4),t97.*t155.*(-t98.*t100.*t156.*ve_d1.*xk1_4+f0v1.*n2k1_2.*t6.*t8.*t80.*ve_v1.*xk1_2+f0d2.*t80.*t85.*t107.*t148.*ve_d2.*xk2_4),t147.*t167.*(-t85.*t107.*t148.*ve_d2.*xk2_4+f0v2.*n2k2_2.*t17.*t19.*t81.*ve_v2.*xk2_2+f0d3.*t81.*t88.*t105.*t150.*ve_d3.*xk3_4),t149.*t169.*(-t88.*t105.*t150.*ve_d3.*xk3_4+f0v3.*n2k3_2.*t25.*t27.*t82.*ve_v3.*xk3_2+f0d4.*n2k4_2.*t35.*t82.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),t151.*t171.*(n2k4_2.*t33.*t35.*ve_v4.*xk4_2-n2k4_2.*t35.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn3_4.^2.*(c33.*t11-c23.*xn3_3)-c33.*t4,t5.*t6.*(n2k1_3-n2k1_3.*t8.*ve_v1),t16.*t17.*(n2k2_3-n2k2_3.*t19.*ve_v2),t24.*t25.*(n2k3_3-n2k3_3.*t27.*ve_v3),t32.*t33.*(n2k4_3-n2k4_3.*t35.*ve_v4),t5.*t40.*(n2k1_3.*t41-n2k1_3.*t6.*t8.*ve_v1.*xk1_2),t16.*t50.*(n2k2_3.*t51-n2k2_3.*t17.*t19.*ve_v2.*xk2_2),t24.*t60.*(n2k3_3.*t61-n2k3_3.*t25.*t27.*ve_v3.*xk3_2),t32.*t70.*(n2k4_3.*t71-n2k4_3.*t33.*t35.*ve_v4.*xk4_2),t97.*t98.*(t109+t116-t100.*t157.*ve_d1),t147.*t148.*(t110+t114-t85.*t115.*ve_d2),t149.*t150.*(t111+t112-t88.*t113.*ve_d3),t151.*t152.*(f0v4.*n2k4_3.*t35.*t153.*ve_v4-n2k4_3.*t35.*t91.*ve_d4.*ve_v4),t97.*t155.*(-t98.*t100.*t157.*ve_d1.*xk1_4+f0v1.*n2k1_3.*t6.*t8.*t80.*ve_v1.*xk1_2+f0d2.*t80.*t85.*t115.*t148.*ve_d2.*xk2_4),t147.*t167.*(-t85.*t115.*t148.*ve_d2.*xk2_4+f0v2.*n2k2_3.*t17.*t19.*t81.*ve_v2.*xk2_2+f0d3.*t81.*t88.*t113.*t150.*ve_d3.*xk3_4),t149.*t169.*(-t88.*t113.*t150.*ve_d3.*xk3_4+f0v3.*n2k3_3.*t25.*t27.*t82.*ve_v3.*xk3_2+f0d4.*n2k4_3.*t35.*t82.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),t151.*t171.*(n2k4_3.*t33.*t35.*ve_v4.*xk4_2-n2k4_3.*t35.*t91.*t152.*ve_d4.*ve_v4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5.*t42.*(t12+t13+t14-t8.*t47+1.0)-t6.*t8.*t15.*t49,0.0,0.0,0.0,t5.*t40.*(t8.*t42.*t47.*xk1_2-t6.*t8.*t15.*t49.*tt_v1.*xk1_2),0.0,0.0,0.0,t97.*t98.*(f0v1.*t8.*t15.*t49.*t80.*tt_v1-f0v1.*t8.*t15.*t49.*t80.*t100.*tt_v1.*ve_d1),0.0,0.0,0.0,-t97.*t155.*(f0v1.*t8.*t42.*t47.*t80.*xk1_2-f0v1.*t6.*t8.*t15.*t49.*t80.*tt_v1.*xk1_2+f0v1.*t8.*t15.*t49.*t80.*t98.*t100.*tt_v1.*ve_d1.*xk1_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t16.*t52.*(t20+t21+t22-t19.*t57+1.0)-t17.*t19.*t23.*t59,0.0,0.0,0.0,t16.*t50.*(t19.*t52.*t57.*xk2_2-t17.*t19.*t23.*t59.*tt_v2.*xk2_2),0.0,0.0,t97.*t98.*(f0v2.*t19.*t23.*t59.*t80.*t85.*tt_v2.*ve_d2-f0v2.*t19.*t23.*t59.*t80.*t85.*t100.*tt_v2.*ve_d1.*ve_d2),t147.*t148.*(f0v2.*t19.*t23.*t59.*t81.*tt_v2-f0v2.*t19.*t23.*t59.*t81.*t85.*tt_v2.*ve_d2),0.0,0.0,t97.*t155.*(f0v2.*t19.*t23.*t59.*t80.*t85.*t148.*tt_v2.*ve_d2.*xk2_4-f0v2.*t19.*t23.*t59.*t80.*t85.*t98.*t100.*tt_v2.*ve_d1.*ve_d2.*xk1_4),-t147.*t167.*(f0v2.*t19.*t52.*t57.*t81.*xk2_2-f0v2.*t17.*t19.*t23.*t59.*t81.*tt_v2.*xk2_2+f0v2.*t19.*t23.*t59.*t81.*t85.*t148.*tt_v2.*ve_d2.*xk2_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t24.*t62.*(t28+t29+t30-t27.*t67+1.0)-t25.*t27.*t31.*t69,0.0,0.0,0.0,t24.*t60.*(t27.*t62.*t67.*xk3_2-t25.*t27.*t31.*t69.*tt_v3.*xk3_2),0.0,t97.*t98.*(f0v3.*t27.*t31.*t69.*t80.*t85.*t88.*tt_v3.*ve_d2.*ve_d3-f0v3.*t27.*t31.*t69.*t80.*t85.*t88.*t100.*tt_v3.*ve_d1.*ve_d2.*ve_d3),t147.*t148.*(f0v3.*t27.*t31.*t69.*t81.*t88.*tt_v3.*ve_d3-f0v3.*t27.*t31.*t69.*t81.*t85.*t88.*tt_v3.*ve_d2.*ve_d3),t149.*t150.*(f0v3.*t27.*t31.*t69.*t82.*tt_v3-f0v3.*t27.*t31.*t69.*t82.*t88.*tt_v3.*ve_d3),0.0,t97.*t155.*(f0v3.*t27.*t31.*t69.*t80.*t85.*t88.*t148.*tt_v3.*ve_d2.*ve_d3.*xk2_4-f0v3.*t27.*t31.*t69.*t80.*t85.*t88.*t98.*t100.*tt_v3.*ve_d1.*ve_d2.*ve_d3.*xk1_4),t147.*t167.*(f0v3.*t27.*t31.*t69.*t81.*t88.*t150.*tt_v3.*ve_d3.*xk3_4-f0v3.*t27.*t31.*t69.*t81.*t85.*t88.*t148.*tt_v3.*ve_d2.*ve_d3.*xk2_4),-t149.*t169.*(f0v3.*t27.*t62.*t67.*t82.*xk3_2-f0v3.*t25.*t27.*t31.*t69.*t82.*tt_v3.*xk3_2+f0v3.*t27.*t31.*t69.*t82.*t88.*t150.*tt_v3.*ve_d3.*xk3_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t32.*t72.*(t36+t37+t38-t35.*t77+1.0)-t33.*t35.*t39.*t79,0.0,0.0,0.0,t32.*t70.*(t172-t33.*t35.*t39.*t79.*tt_v4.*xk4_2),t97.*t98.*(f0d4.*t35.*t39.*t79.*t80.*t85.*t88.*t91.*tt_v4.*ve_d2.*ve_d3.*ve_d4-f0d4.*t35.*t39.*t79.*t80.*t85.*t88.*t91.*t100.*tt_v4.*ve_d1.*ve_d2.*ve_d3.*ve_d4),t147.*t148.*(f0d4.*t35.*t39.*t79.*t81.*t88.*t91.*tt_v4.*ve_d3.*ve_d4-f0d4.*t35.*t39.*t79.*t81.*t85.*t88.*t91.*tt_v4.*ve_d2.*ve_d3.*ve_d4),t149.*t150.*(f0d4.*t35.*t39.*t79.*t82.*t91.*tt_v4.*ve_d4-f0d4.*t35.*t39.*t79.*t82.*t88.*t91.*tt_v4.*ve_d3.*ve_d4),t151.*t152.*(f0v4.*t35.*t39.*t79.*t153.*tt_v4-t35.*t39.*t79.*t91.*tt_v4.*ve_d4),t97.*t155.*(f0d4.*t35.*t39.*t79.*t80.*t85.*t88.*t91.*t148.*tt_v4.*ve_d2.*ve_d3.*ve_d4.*xk2_4-f0d4.*t35.*t39.*t79.*t80.*t85.*t88.*t91.*t98.*t100.*tt_v4.*ve_d1.*ve_d2.*ve_d3.*ve_d4.*xk1_4),t147.*t167.*(f0d4.*t35.*t39.*t79.*t81.*t88.*t91.*t150.*tt_v4.*ve_d3.*ve_d4.*xk3_4-f0d4.*t35.*t39.*t79.*t81.*t85.*t88.*t91.*t148.*tt_v4.*ve_d2.*ve_d3.*ve_d4.*xk2_4),t149.*t169.*(f0d4.*t35.*t39.*t79.*t82.*t91.*t152.*tt_v4.*ve_d4.*xk4_4-f0d4.*t35.*t39.*t79.*t82.*t88.*t91.*t150.*tt_v4.*ve_d3.*ve_d4.*xk3_4),-t151.*t171.*(t172-t33.*t35.*t39.*t79.*tt_v4.*xk4_2+t35.*t39.*t79.*t91.*t152.*tt_v4.*ve_d4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5.*1.0./xk1_2.^2.*(t41.*(nr1+t12+t13+t14)-t6.*t8.*t47.*xk1_2)-t5.*t6.*t8.*t40.*t47,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v1.*t6.*t8.*t47.*t80.*t97.*t155,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t16.*1.0./xk2_2.^2.*(t51.*(nr2+t20+t21+t22)-t17.*t19.*t57.*xk2_2)-t16.*t17.*t19.*t50.*t57,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v2.*t17.*t19.*t57.*t81.*t147.*t167,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t24.*1.0./xk3_2.^2.*(t61.*(nr3+t28+t29+t30)-t25.*t27.*t67.*xk3_2)-t24.*t25.*t27.*t60.*t67,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v3.*t25.*t27.*t67.*t82.*t149.*t169,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t32.*1.0./xk4_2.^2.*(t71.*(nr4+t36+t37+t38)-t33.*t35.*t77.*xk4_2)-t32.*t33.*t35.*t70.*t77,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t33.*t35.*t77.*t151.*t171,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t97.*t158.*(t117+t139-t100.*t163)-t98.*t100.*t140.*t165,0.0,0.0,0.0,t97.*t155.*(t100.*t158.*t163.*xk1_4-t98.*t100.*t140.*t165.*tt_d1.*xk1_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t97.*t98.*(f0d2.*t80.*t85.*t135.*t142.*tt_d2-f0d2.*t80.*t85.*t100.*t135.*t142.*tt_d2.*ve_d1),-t147.*t166.*(t118+t132-t85.*t138)-t85.*t135.*t142.*t148,0.0,0.0,-t97.*t155.*(f0d2.*t80.*t85.*t138.*t166.*xk2_4-f0d2.*t80.*t85.*t135.*t142.*t148.*tt_d2.*xk2_4+f0d2.*t80.*t85.*t98.*t100.*t135.*t142.*tt_d2.*ve_d1.*xk1_4),t147.*t167.*(t85.*t138.*t166.*xk2_4-t85.*t135.*t142.*t148.*tt_d2.*xk2_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t97.*t98.*(f0d3.*t80.*t85.*t88.*t119.*t144.*tt_d3.*ve_d2-f0d3.*t80.*t85.*t88.*t100.*t119.*t144.*tt_d3.*ve_d1.*ve_d2),t147.*t148.*(f0d3.*t81.*t88.*t119.*t144.*tt_d3-f0d3.*t81.*t85.*t88.*t119.*t144.*tt_d3.*ve_d2),-t149.*t168.*(t122+t128-t88.*t131)-t88.*t119.*t144.*t150,0.0,t97.*t155.*(f0d3.*t80.*t85.*t88.*t119.*t144.*t148.*tt_d3.*ve_d2.*xk2_4-f0d3.*t80.*t85.*t88.*t98.*t100.*t119.*t144.*tt_d3.*ve_d1.*ve_d2.*xk1_4),-t147.*t167.*(f0d3.*t81.*t88.*t131.*t168.*xk3_4-f0d3.*t81.*t88.*t119.*t144.*t150.*tt_d3.*xk3_4+f0d3.*t81.*t85.*t88.*t119.*t144.*t148.*tt_d3.*ve_d2.*xk2_4),t149.*t169.*(t88.*t131.*t168.*xk3_4-t88.*t119.*t144.*t150.*tt_d3.*xk3_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t97.*t98.*(f0d4.*t80.*t85.*t88.*t91.*t123.*t146.*tt_d4.*ve_d2.*ve_d3-f0d4.*t80.*t85.*t88.*t91.*t100.*t123.*t146.*tt_d4.*ve_d1.*ve_d2.*ve_d3),t147.*t148.*(f0d4.*t81.*t88.*t91.*t123.*t146.*tt_d4.*ve_d3-f0d4.*t81.*t85.*t88.*t91.*t123.*t146.*tt_d4.*ve_d2.*ve_d3),t149.*t150.*(f0d4.*t82.*t91.*t123.*t146.*tt_d4-f0d4.*t82.*t88.*t91.*t123.*t146.*tt_d4.*ve_d3),t151.*t170.*(t91.*t127-f0v4.*t35.*t77.*t153)-t91.*t123.*t146.*t152,t97.*t155.*(f0d4.*t80.*t85.*t88.*t91.*t123.*t146.*t148.*tt_d4.*ve_d2.*ve_d3.*xk2_4-f0d4.*t80.*t85.*t88.*t91.*t98.*t100.*t123.*t146.*tt_d4.*ve_d1.*ve_d2.*ve_d3.*xk1_4),t147.*t167.*(f0d4.*t81.*t88.*t91.*t123.*t146.*t150.*tt_d4.*ve_d3.*xk3_4-f0d4.*t81.*t85.*t88.*t91.*t123.*t146.*t148.*tt_d4.*ve_d2.*ve_d3.*xk2_4),-t149.*t169.*(f0d4.*t82.*t91.*t127.*t170.*xk4_4-f0d4.*t82.*t91.*t123.*t146.*t152.*tt_d4.*xk4_4+f0d4.*t82.*t88.*t91.*t123.*t146.*t150.*tt_d4.*ve_d3.*xk3_4),t151.*t171.*(t91.*t127.*t170.*xk4_4-t91.*t123.*t146.*t152.*tt_d4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t97.*1.0./xk1_4.^2.*(-t98.*t100.*t163.*xk1_4+f0d2.*t80.*t85.*t138.*t148.*xk2_4+f0v1.*t6.*t8.*t47.*t80.*xk1_2)-t97.*t98.*t100.*t155.*t163,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d2.*t80.*t85.*t97.*t138.*t148.*t155,-t147.*1.0./xk2_4.^2.*(-t85.*t138.*t148.*xk2_4+f0d3.*t81.*t88.*t131.*t150.*xk3_4+f0v2.*t17.*t19.*t57.*t81.*xk2_2)-t85.*t138.*t147.*t148.*t167,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d3.*t81.*t88.*t131.*t147.*t150.*t167,-t149.*1.0./xk3_4.^2.*(-t88.*t131.*t150.*xk3_4+f0d4.*t82.*t91.*t127.*t152.*xk4_4+f0v3.*t25.*t27.*t67.*t82.*xk3_2)-t88.*t131.*t149.*t150.*t169,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d4.*t82.*t91.*t127.*t149.*t152.*t169,-t151.*1.0./xk4_4.^2.*(t33.*t35.*t77.*xk4_2-t91.*t127.*t152.*xk4_4)-t91.*t127.*t151.*t152.*t171],[28,28]);
