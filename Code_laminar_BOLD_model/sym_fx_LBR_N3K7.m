function f_s = sym_fx_LBR_N3K7(xn1_1,xn2_1,xn3_1,xn1_2,xn2_2,xn3_2,xn1_3,xn2_3,xn3_3,xn1_4,xn2_4,xn3_4,xk1_1,xk2_1,xk3_1,xk4_1,xk5_1,xk6_1,xk7_1,xk1_2,xk2_2,xk3_2,xk4_2,xk5_2,xk6_2,xk7_2,xk1_3,xk2_3,xk3_3,xk4_3,xk5_3,xk6_3,xk7_3,xk1_4,xk2_4,xk3_4,xk4_4,xk5_4,xk6_4,xk7_4,a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3,mu1,mu2,mu3,cu1,cu2,cu3,lam1,lam2,lam3,c11,c12,c13,c21,c22,c23,c31,c32,c33,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k5_1,n2k6_1,n2k7_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k5_2,n2k6_2,n2k7_2,n2k1_3,n2k2_3,n2k3_3,n2k4_3,n2k5_3,n2k6_3,n2k7_3,al_v1,al_v2,al_v3,al_v4,al_v5,al_v6,al_v7,al_d1,al_d2,al_d3,al_d4,al_d5,al_d6,al_d7,tt_v1,tt_v2,tt_v3,tt_v4,tt_v5,tt_v6,tt_v7,tt_d1,tt_d2,tt_d3,tt_d4,tt_d5,tt_d6,tt_d7,nr1,nr2,nr3,nr4,nr5,nr6,nr7,ve_v1,ve_v2,ve_v3,ve_v4,ve_v5,ve_v6,ve_v7,ve_d1,ve_d2,ve_d3,ve_d4,ve_d5,ve_d6,ve_d7,f0v1,f0v2,f0v3,f0v4,f0v5,f0v6,f0v7,f0d1,f0d2,f0d3,f0d4,f0d5,f0d6,f0d7)
%SYM_FX_LBR_N3K7
%    F_S = SYM_FX_LBR_N3K7(XN1_1,XN2_1,XN3_1,XN1_2,XN2_2,XN3_2,XN1_3,XN2_3,XN3_3,XN1_4,XN2_4,XN3_4,XK1_1,XK2_1,XK3_1,XK4_1,XK5_1,XK6_1,XK7_1,XK1_2,XK2_2,XK3_2,XK4_2,XK5_2,XK6_2,XK7_2,XK1_3,XK2_3,XK3_3,XK4_3,XK5_3,XK6_3,XK7_3,XK1_4,XK2_4,XK3_4,XK4_4,XK5_4,XK6_4,XK7_4,A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,MU1,MU2,MU3,CU1,CU2,CU3,LAM1,LAM2,LAM3,C11,C12,C13,C21,C22,C23,C31,C32,C33,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K5_1,N2K6_1,N2K7_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K5_2,N2K6_2,N2K7_2,N2K1_3,N2K2_3,N2K3_3,N2K4_3,N2K5_3,N2K6_3,N2K7_3,AL_V1,AL_V2,AL_V3,AL_V4,AL_V5,AL_V6,AL_V7,AL_D1,AL_D2,AL_D3,AL_D4,AL_D5,AL_D6,AL_D7,TT_V1,TT_V2,TT_V3,TT_V4,TT_V5,TT_V6,TT_V7,TT_D1,TT_D2,TT_D3,TT_D4,TT_D5,TT_D6,TT_D7,NR1,NR2,NR3,NR4,NR5,NR6,NR7,VE_V1,VE_V2,VE_V3,VE_V4,VE_V5,VE_V6,VE_V7,VE_D1,VE_D2,VE_D3,VE_D4,VE_D5,VE_D6,VE_D7,F0V1,F0V2,F0V3,F0V4,F0V5,F0V6,F0V7,F0D1,F0D2,F0D3,F0D4,F0D5,F0D6,F0D7)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    10-Mar-2020 11:08:38

t2 = xn1_4-1.0;
t3 = xn2_4-1.0;
t4 = xn3_4-1.0;
t5 = n2k1_1.*t2;
t6 = n2k1_2.*t3;
t7 = n2k1_3.*t4;
t8 = n2k2_1.*t2;
t9 = n2k2_2.*t3;
t10 = n2k2_3.*t4;
t11 = n2k3_1.*t2;
t12 = n2k3_2.*t3;
t13 = n2k3_3.*t4;
t14 = n2k4_1.*t2;
t15 = n2k4_2.*t3;
t16 = n2k4_3.*t4;
t17 = n2k5_1.*t2;
t18 = n2k5_2.*t3;
t19 = n2k5_3.*t4;
t20 = n2k6_1.*t2;
t21 = n2k6_2.*t3;
t22 = n2k6_3.*t4;
t23 = n2k7_1.*t2;
t24 = n2k7_2.*t3;
t25 = n2k7_3.*t4;
t26 = 1.0./tt_v1;
t27 = 1.0./xk1_1;
t28 = tt_v1+ve_v1;
t29 = 1.0./t28;
t30 = t5+t6+t7+1.0;
t31 = t30.*ve_v1;
t32 = 1.0./al_v1;
t33 = xk1_1.^t32;
t34 = t33.*tt_v1;
t35 = t31+t34;
t36 = 1.0./tt_v2;
t37 = 1.0./xk2_1;
t38 = tt_v2+ve_v2;
t39 = 1.0./t38;
t40 = t8+t9+t10+1.0;
t41 = t40.*ve_v2;
t42 = 1.0./al_v2;
t43 = xk2_1.^t42;
t44 = t43.*tt_v2;
t45 = t41+t44;
t46 = 1.0./tt_v3;
t47 = 1.0./xk3_1;
t48 = tt_v3+ve_v3;
t49 = 1.0./t48;
t50 = t11+t12+t13+1.0;
t51 = t50.*ve_v3;
t52 = 1.0./al_v3;
t53 = xk3_1.^t52;
t54 = t53.*tt_v3;
t55 = t51+t54;
t56 = 1.0./tt_v4;
t57 = 1.0./xk4_1;
t58 = tt_v4+ve_v4;
t59 = 1.0./t58;
t60 = t14+t15+t16+1.0;
t61 = t60.*ve_v4;
t62 = 1.0./al_v4;
t63 = xk4_1.^t62;
t64 = t63.*tt_v4;
t65 = t61+t64;
t66 = 1.0./tt_v5;
t67 = 1.0./xk5_1;
t68 = tt_v5+ve_v5;
t69 = 1.0./t68;
t70 = t17+t18+t19+1.0;
t71 = t70.*ve_v5;
t72 = 1.0./al_v5;
t73 = xk5_1.^t72;
t74 = t73.*tt_v5;
t75 = t71+t74;
t76 = 1.0./tt_v6;
t77 = 1.0./xk6_1;
t78 = tt_v6+ve_v6;
t79 = 1.0./t78;
t80 = t20+t21+t22+1.0;
t81 = t80.*ve_v6;
t82 = 1.0./al_v6;
t83 = xk6_1.^t82;
t84 = t83.*tt_v6;
t85 = t81+t84;
t86 = 1.0./tt_v7;
t87 = 1.0./xk7_1;
t88 = tt_v7+ve_v7;
t89 = 1.0./t88;
t90 = t23+t24+t25+1.0;
t91 = t90.*ve_v7;
t92 = 1.0./al_v7;
t93 = xk7_1.^t92;
t94 = t93.*tt_v7;
t95 = t91+t94;
t96 = 1.0./f0d1;
t97 = 1.0./f0d2;
t98 = 1.0./f0d3;
t99 = 1.0./f0d4;
t100 = 1.0./f0d5;
t101 = 1.0./f0d6;
t102 = f0v1.*t29.*t35.*t96;
t103 = f0v2.*t39.*t45.*t97;
t104 = tt_d3+ve_d3;
t105 = 1.0./t104;
t106 = 1.0./al_d3;
t107 = xk3_3.^t106;
t108 = t107.*tt_d3;
t109 = f0v3.*t49.*t55.*t98;
t110 = tt_d4+ve_d4;
t111 = 1.0./t110;
t112 = 1.0./al_d4;
t113 = xk4_3.^t112;
t114 = t113.*tt_d4;
t115 = f0v4.*t59.*t65.*t99;
t116 = tt_d5+ve_d5;
t117 = 1.0./t116;
t118 = f0v5.*t69.*t75.*t100;
t119 = tt_d6+ve_d6;
t120 = 1.0./t119;
t121 = 1.0./al_d6;
t122 = xk6_3.^t121;
t123 = t122.*tt_d6;
t124 = f0v6.*t79.*t85.*t101;
t125 = 1.0./al_d7;
t126 = xk7_3.^t125;
t127 = t126.*tt_d7;
t128 = t89.*t95.*ve_d7;
t129 = t127+t128;
t130 = tt_d7+ve_d7;
t131 = 1.0./t130;
t132 = f0d7.*t101.*t129.*t131;
t133 = t124+t132;
t134 = t133.*ve_d6;
t135 = t123+t134;
t136 = f0d6.*t100.*t120.*t135;
t137 = t118+t136;
t138 = t137.*ve_d5;
t139 = 1.0./al_d5;
t140 = xk5_3.^t139;
t141 = t140.*tt_d5;
t142 = t138+t141;
t143 = f0d5.*t99.*t117.*t142;
t144 = t115+t143;
t145 = t144.*ve_d4;
t146 = t114+t145;
t147 = f0d4.*t98.*t111.*t146;
t148 = t109+t147;
t149 = t148.*ve_d3;
t150 = t108+t149;
t151 = f0d3.*t97.*t105.*t150;
t152 = t103+t151;
t153 = t152.*ve_d2;
t154 = 1.0./al_d2;
t155 = xk2_3.^t154;
t156 = t155.*tt_d2;
t157 = t153+t156;
t158 = tt_d2+ve_d2;
t159 = 1.0./t158;
t160 = f0d2.*t96.*t157.*t159;
t161 = 1.0./tt_d1;
t162 = 1.0./xk1_3;
t163 = tt_d1+ve_d1;
t164 = 1.0./t163;
t165 = 1.0./al_d1;
t166 = xk1_3.^t165;
t167 = t166.*tt_d1;
t168 = t102+t160;
t169 = t168.*ve_d1;
t170 = t167+t169;
t171 = 1.0./xk2_3;
t172 = 1.0./tt_d2;
t173 = 1.0./xk3_3;
t174 = 1.0./tt_d3;
t175 = 1.0./xk4_3;
t176 = 1.0./tt_d4;
t177 = 1.0./xk5_3;
t178 = 1.0./tt_d5;
t179 = 1.0./xk6_3;
t180 = 1.0./tt_d6;
t181 = 1.0./xk7_3;
t182 = 1.0./tt_d7;
f_s = [cu1+a1_1.*xn1_1+a1_2.*xn2_1+a1_3.*xn3_1-mu1.*xn1_2;cu2+a2_1.*xn1_1+a2_2.*xn2_1+a2_3.*xn3_1-mu2.*xn2_2;cu3+a3_1.*xn1_1+a3_2.*xn2_1+a3_3.*xn3_1-mu3.*xn3_2;lam1.*(xn1_1-xn1_2);lam2.*(xn2_1-xn2_2);lam3.*(xn3_1-xn3_2);xn1_1-c11.*xn1_3;xn2_1-c12.*xn2_3;xn3_1-c13.*xn3_3;-(c31.*t2-c21.*xn1_3)./xn1_4;-(c32.*t3-c22.*xn2_3)./xn2_4;-(c33.*t4-c23.*xn3_3)./xn3_4;t26.*t27.*(t5+t6+t7-t29.*t35+1.0);t36.*t37.*(t8+t9+t10-t39.*t45+1.0);t46.*t47.*(t11+t12+t13-t49.*t55+1.0);t56.*t57.*(t14+t15+t16-t59.*t65+1.0);t66.*t67.*(t17+t18+t19-t69.*t75+1.0);t76.*t77.*(t20+t21+t22-t79.*t85+1.0);t86.*t87.*(t23+t24+t25-t89.*t95+1.0);(t26.*((nr1+t5+t6+t7)./nr1-t27.*t29.*t35.*xk1_2))./xk1_2;(t36.*((nr2+t8+t9+t10)./nr2-t37.*t39.*t45.*xk2_2))./xk2_2;(t46.*((nr3+t11+t12+t13)./nr3-t47.*t49.*t55.*xk3_2))./xk3_2;(t56.*((nr4+t14+t15+t16)./nr4-t57.*t59.*t65.*xk4_2))./xk4_2;(t66.*((nr5+t17+t18+t19)./nr5-t67.*t69.*t75.*xk5_2))./xk5_2;(t76.*((nr6+t20+t21+t22)./nr6-t77.*t79.*t85.*xk6_2))./xk6_2;(t86.*((nr7+t23+t24+t25)./nr7-t87.*t89.*t95.*xk7_2))./xk7_2;t161.*t162.*(t102+t160-t164.*t170);t171.*t172.*(t103+t151-t157.*t159);t173.*t174.*(t109+t147-t105.*t150);t175.*t176.*(t115+t143-t111.*t146);t177.*t178.*(t118+t136-t117.*t142);t179.*t180.*(t124+t132-t120.*t135);-t181.*t182.*(t129.*t131-(f0v7.*t89.*t95)./f0d7);(t161.*(-t162.*t164.*t170.*xk1_4+f0d2.*t96.*t157.*t159.*t171.*xk2_4+f0v1.*t27.*t29.*t35.*t96.*xk1_2))./xk1_4;(t172.*(-t157.*t159.*t171.*xk2_4+f0d3.*t97.*t105.*t150.*t173.*xk3_4+f0v2.*t37.*t39.*t45.*t97.*xk2_2))./xk2_4;(t174.*(-t105.*t150.*t173.*xk3_4+f0d4.*t98.*t111.*t146.*t175.*xk4_4+f0v3.*t47.*t49.*t55.*t98.*xk3_2))./xk3_4;(t176.*(-t111.*t146.*t175.*xk4_4+f0d5.*t99.*t117.*t142.*t177.*xk5_4+f0v4.*t57.*t59.*t65.*t99.*xk4_2))./xk4_4;(t178.*(-t117.*t142.*t177.*xk5_4+f0d6.*t100.*t120.*t135.*t179.*xk6_4+f0v5.*t67.*t69.*t75.*t100.*xk5_2))./xk5_4;(t180.*(-t120.*t135.*t179.*xk6_4+f0d7.*t101.*t129.*t131.*t181.*xk7_4+f0v6.*t77.*t79.*t85.*t101.*xk6_2))./xk6_4;(t182.*(t87.*t89.*t95.*xk7_2-t129.*t131.*t181.*xk7_4))./xk7_4];
