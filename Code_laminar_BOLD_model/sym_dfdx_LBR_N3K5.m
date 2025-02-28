function dfdx_s = sym_dfdx_LBR_N3K5(xn1_1,xn2_1,xn3_1,xn1_2,xn2_2,xn3_2,xn1_3,xn2_3,xn3_3,xn1_4,xn2_4,xn3_4,xk1_1,xk2_1,xk3_1,xk4_1,xk5_1,xk1_2,xk2_2,xk3_2,xk4_2,xk5_2,xk1_3,xk2_3,xk3_3,xk4_3,xk5_3,xk1_4,xk2_4,xk3_4,xk4_4,xk5_4,a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3,mu1,mu2,mu3,cu1,cu2,cu3,lam1,lam2,lam3,c11,c12,c13,c21,c22,c23,c31,c32,c33,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k5_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k5_2,n2k1_3,n2k2_3,n2k3_3,n2k4_3,n2k5_3,al_v1,al_v2,al_v3,al_v4,al_v5,al_d1,al_d2,al_d3,al_d4,al_d5,tt_v1,tt_v2,tt_v3,tt_v4,tt_v5,tt_d1,tt_d2,tt_d3,tt_d4,tt_d5,nr1,nr2,nr3,nr4,nr5,ve_v1,ve_v2,ve_v3,ve_v4,ve_v5,ve_d1,ve_d2,ve_d3,ve_d4,ve_d5,f0v1,f0v2,f0v3,f0v4,f0v5,f0d1,f0d2,f0d3,f0d4,f0d5)
%SYM_DFDX_LBR_N3K5
%    DFDX_S = SYM_DFDX_LBR_N3K5(XN1_1,XN2_1,XN3_1,XN1_2,XN2_2,XN3_2,XN1_3,XN2_3,XN3_3,XN1_4,XN2_4,XN3_4,XK1_1,XK2_1,XK3_1,XK4_1,XK5_1,XK1_2,XK2_2,XK3_2,XK4_2,XK5_2,XK1_3,XK2_3,XK3_3,XK4_3,XK5_3,XK1_4,XK2_4,XK3_4,XK4_4,XK5_4,A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,MU1,MU2,MU3,CU1,CU2,CU3,LAM1,LAM2,LAM3,C11,C12,C13,C21,C22,C23,C31,C32,C33,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K5_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K5_2,N2K1_3,N2K2_3,N2K3_3,N2K4_3,N2K5_3,AL_V1,AL_V2,AL_V3,AL_V4,AL_V5,AL_D1,AL_D2,AL_D3,AL_D4,AL_D5,TT_V1,TT_V2,TT_V3,TT_V4,TT_V5,TT_D1,TT_D2,TT_D3,TT_D4,TT_D5,NR1,NR2,NR3,NR4,NR5,VE_V1,VE_V2,VE_V3,VE_V4,VE_V5,VE_D1,VE_D2,VE_D3,VE_D4,VE_D5,F0V1,F0V2,F0V3,F0V4,F0V5,F0D1,F0D2,F0D3,F0D4,F0D5)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    10-Mar-2020 11:04:39

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
t40 = 1.0./tt_v5;
t41 = 1.0./xk5_1;
t42 = tt_v5+ve_v5;
t43 = 1.0./t42;
t44 = n2k5_1.*t9;
t45 = n2k5_2.*t10;
t46 = n2k5_3.*t11;
t47 = 1.0./al_v5;
t48 = 1.0./xk1_2;
t49 = 1.0./nr1;
t50 = 1.0./xk1_1.^2;
t51 = t12+t13+t14+1.0;
t52 = t51.*ve_v1;
t53 = xk1_1.^t15;
t54 = t53.*tt_v1;
t55 = t52+t54;
t56 = t15-1.0;
t57 = xk1_1.^t56;
t58 = 1.0./xk2_2;
t59 = 1.0./nr2;
t60 = 1.0./xk2_1.^2;
t61 = t20+t21+t22+1.0;
t62 = t61.*ve_v2;
t63 = xk2_1.^t23;
t64 = t63.*tt_v2;
t65 = t62+t64;
t66 = t23-1.0;
t67 = xk2_1.^t66;
t68 = 1.0./xk3_2;
t69 = 1.0./nr3;
t70 = 1.0./xk3_1.^2;
t71 = t28+t29+t30+1.0;
t72 = t71.*ve_v3;
t73 = xk3_1.^t31;
t74 = t73.*tt_v3;
t75 = t72+t74;
t76 = t31-1.0;
t77 = xk3_1.^t76;
t78 = 1.0./xk4_2;
t79 = 1.0./nr4;
t80 = 1.0./xk4_1.^2;
t81 = t36+t37+t38+1.0;
t82 = t81.*ve_v4;
t83 = xk4_1.^t39;
t84 = t83.*tt_v4;
t85 = t82+t84;
t86 = t39-1.0;
t87 = xk4_1.^t86;
t88 = 1.0./xk5_2;
t89 = 1.0./nr5;
t90 = 1.0./xk5_1.^2;
t91 = t44+t45+t46+1.0;
t92 = t91.*ve_v5;
t93 = xk5_1.^t47;
t94 = t93.*tt_v5;
t95 = t92+t94;
t96 = t47-1.0;
t97 = xk5_1.^t96;
t98 = 1.0./f0d1;
t99 = 1.0./f0d2;
t100 = 1.0./f0d3;
t101 = 1.0./f0d4;
t102 = f0v1.*n2k1_1.*t8.*t98.*ve_v1;
t103 = tt_d2+ve_d2;
t104 = 1.0./t103;
t105 = f0v2.*n2k2_1.*t19.*t99.*ve_v2;
t106 = tt_d3+ve_d3;
t107 = 1.0./t106;
t108 = f0v3.*n2k3_1.*t27.*t100.*ve_v3;
t109 = tt_d4+ve_d4;
t110 = 1.0./t109;
t111 = f0v4.*n2k4_1.*t35.*t101.*ve_v4;
t112 = tt_d5+ve_d5;
t113 = 1.0./t112;
t114 = f0d5.*n2k5_1.*t43.*t101.*t113.*ve_d5.*ve_v5;
t115 = t111+t114;
t116 = f0d4.*t100.*t110.*t115.*ve_d4;
t117 = t108+t116;
t118 = f0d3.*t99.*t107.*t117.*ve_d3;
t119 = t105+t118;
t120 = f0d2.*t98.*t104.*t119.*ve_d2;
t121 = 1.0./tt_d1;
t122 = 1.0./xk1_3;
t123 = tt_d1+ve_d1;
t124 = 1.0./t123;
t125 = f0v1.*n2k1_2.*t8.*t98.*ve_v1;
t126 = f0v2.*n2k2_2.*t19.*t99.*ve_v2;
t127 = f0v3.*n2k3_2.*t27.*t100.*ve_v3;
t128 = f0v4.*n2k4_2.*t35.*t101.*ve_v4;
t129 = f0d5.*n2k5_2.*t43.*t101.*t113.*ve_d5.*ve_v5;
t130 = t128+t129;
t131 = f0d4.*t100.*t110.*t130.*ve_d4;
t132 = t127+t131;
t133 = f0d3.*t99.*t107.*t132.*ve_d3;
t134 = t126+t133;
t135 = f0d2.*t98.*t104.*t134.*ve_d2;
t136 = f0v1.*n2k1_3.*t8.*t98.*ve_v1;
t137 = f0v2.*n2k2_3.*t19.*t99.*ve_v2;
t138 = f0v3.*n2k3_3.*t27.*t100.*ve_v3;
t139 = f0v4.*n2k4_3.*t35.*t101.*ve_v4;
t140 = f0d5.*n2k5_3.*t43.*t101.*t113.*ve_d5.*ve_v5;
t141 = t139+t140;
t142 = f0d4.*t100.*t110.*t141.*ve_d4;
t143 = t138+t142;
t144 = f0d3.*t99.*t107.*t143.*ve_d3;
t145 = t137+t144;
t146 = f0d2.*t98.*t104.*t145.*ve_d2;
t147 = f0v1.*t8.*t55.*t98;
t148 = 1.0./al_d2;
t149 = xk2_3.^t148;
t150 = t149.*tt_d2;
t151 = f0v2.*t19.*t65.*t99;
t152 = f0v3.*t27.*t75.*t100;
t153 = 1.0./al_d4;
t154 = xk4_3.^t153;
t155 = t154.*tt_d4;
t156 = f0v4.*t35.*t85.*t101;
t157 = 1.0./al_d5;
t158 = xk5_3.^t157;
t159 = t158.*tt_d5;
t160 = t43.*t95.*ve_d5;
t161 = t159+t160;
t162 = f0d5.*t101.*t113.*t161;
t163 = t156+t162;
t164 = t163.*ve_d4;
t165 = t155+t164;
t166 = f0d4.*t100.*t110.*t165;
t167 = t152+t166;
t168 = t167.*ve_d3;
t169 = 1.0./al_d3;
t170 = xk3_3.^t169;
t171 = t170.*tt_d3;
t172 = t168+t171;
t173 = f0d3.*t99.*t107.*t172;
t174 = t151+t173;
t175 = t174.*ve_d2;
t176 = t150+t175;
t177 = f0d2.*t98.*t104.*t176;
t178 = 1.0./al_d1;
t179 = t148-1.0;
t180 = xk2_3.^t179;
t181 = t169-1.0;
t182 = xk3_3.^t181;
t183 = t153-1.0;
t184 = xk4_3.^t183;
t185 = t157-1.0;
t186 = xk5_3.^t185;
t187 = 1.0./tt_d2;
t188 = 1.0./xk2_3;
t189 = 1.0./tt_d3;
t190 = 1.0./xk3_3;
t191 = 1.0./tt_d4;
t192 = 1.0./xk4_3;
t193 = 1.0./tt_d5;
t194 = 1.0./xk5_3;
t195 = 1.0./f0d5;
t196 = t102+t120;
t197 = 1.0./xk1_4;
t198 = t125+t135;
t199 = t136+t146;
t200 = 1.0./xk1_3.^2;
t201 = xk1_3.^t178;
t202 = t201.*tt_d1;
t203 = t147+t177;
t204 = t203.*ve_d1;
t205 = t202+t204;
t206 = t178-1.0;
t207 = xk1_3.^t206;
t208 = 1.0./xk2_3.^2;
t209 = 1.0./xk2_4;
t210 = 1.0./xk3_3.^2;
t211 = 1.0./xk3_4;
t212 = 1.0./xk4_3.^2;
t213 = 1.0./xk4_4;
t214 = 1.0./xk5_3.^2;
t215 = 1.0./xk5_4;
t216 = t43.*t90.*t95.*xk5_2;
dfdx_s = reshape([a1_1,a2_1,a3_1,lam1,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,a1_2,a2_2,a3_2,0.0,lam2,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,a1_3,a2_3,a3_3,0.0,0.0,lam3,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu1,0.0,0.0,-lam1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu2,0.0,0.0,-lam2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-mu3,0.0,0.0,-lam3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c11,0.0,0.0,c21.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c12,0.0,0.0,c22.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-c13,0.0,0.0,c23.*t4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn1_4.^2.*(c31.*t9-c21.*xn1_3)-c31.*t2,0.0,0.0,t5.*t6.*(n2k1_1-n2k1_1.*t8.*ve_v1),t16.*t17.*(n2k2_1-n2k2_1.*t19.*ve_v2),t24.*t25.*(n2k3_1-n2k3_1.*t27.*ve_v3),t32.*t33.*(n2k4_1-n2k4_1.*t35.*ve_v4),t40.*t41.*(n2k5_1-n2k5_1.*t43.*ve_v5),t5.*t48.*(n2k1_1.*t49-n2k1_1.*t6.*t8.*ve_v1.*xk1_2),t16.*t58.*(n2k2_1.*t59-n2k2_1.*t17.*t19.*ve_v2.*xk2_2),t24.*t68.*(n2k3_1.*t69-n2k3_1.*t25.*t27.*ve_v3.*xk3_2),t32.*t78.*(n2k4_1.*t79-n2k4_1.*t33.*t35.*ve_v4.*xk4_2),t40.*t88.*(n2k5_1.*t89-n2k5_1.*t41.*t43.*ve_v5.*xk5_2),t121.*t122.*(t102+t120-t124.*t196.*ve_d1),t187.*t188.*(t105+t118-t104.*t119.*ve_d2),t189.*t190.*(t108+t116-t107.*t117.*ve_d3),t191.*t192.*(t111+t114-t110.*t115.*ve_d4),t193.*t194.*(f0v5.*n2k5_1.*t43.*t195.*ve_v5-n2k5_1.*t43.*t113.*ve_d5.*ve_v5),t121.*t197.*(-t122.*t124.*t196.*ve_d1.*xk1_4+f0v1.*n2k1_1.*t6.*t8.*t98.*ve_v1.*xk1_2+f0d2.*t98.*t104.*t119.*t188.*ve_d2.*xk2_4),t187.*t209.*(-t104.*t119.*t188.*ve_d2.*xk2_4+f0v2.*n2k2_1.*t17.*t19.*t99.*ve_v2.*xk2_2+f0d3.*t99.*t107.*t117.*t190.*ve_d3.*xk3_4),t189.*t211.*(-t107.*t117.*t190.*ve_d3.*xk3_4+f0v3.*n2k3_1.*t25.*t27.*t100.*ve_v3.*xk3_2+f0d4.*t100.*t110.*t115.*t192.*ve_d4.*xk4_4),t191.*t213.*(-t110.*t115.*t192.*ve_d4.*xk4_4+f0v4.*n2k4_1.*t33.*t35.*t101.*ve_v4.*xk4_2+f0d5.*n2k5_1.*t43.*t101.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),t193.*t215.*(n2k5_1.*t41.*t43.*ve_v5.*xk5_2-n2k5_1.*t43.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn2_4.^2.*(c32.*t10-c22.*xn2_3)-c32.*t3,0.0,t5.*t6.*(n2k1_2-n2k1_2.*t8.*ve_v1),t16.*t17.*(n2k2_2-n2k2_2.*t19.*ve_v2),t24.*t25.*(n2k3_2-n2k3_2.*t27.*ve_v3),t32.*t33.*(n2k4_2-n2k4_2.*t35.*ve_v4),t40.*t41.*(n2k5_2-n2k5_2.*t43.*ve_v5),t5.*t48.*(n2k1_2.*t49-n2k1_2.*t6.*t8.*ve_v1.*xk1_2),t16.*t58.*(n2k2_2.*t59-n2k2_2.*t17.*t19.*ve_v2.*xk2_2),t24.*t68.*(n2k3_2.*t69-n2k3_2.*t25.*t27.*ve_v3.*xk3_2),t32.*t78.*(n2k4_2.*t79-n2k4_2.*t33.*t35.*ve_v4.*xk4_2),t40.*t88.*(n2k5_2.*t89-n2k5_2.*t41.*t43.*ve_v5.*xk5_2),t121.*t122.*(t125+t135-t124.*t198.*ve_d1),t187.*t188.*(t126+t133-t104.*t134.*ve_d2),t189.*t190.*(t127+t131-t107.*t132.*ve_d3),t191.*t192.*(t128+t129-t110.*t130.*ve_d4),t193.*t194.*(f0v5.*n2k5_2.*t43.*t195.*ve_v5-n2k5_2.*t43.*t113.*ve_d5.*ve_v5),t121.*t197.*(-t122.*t124.*t198.*ve_d1.*xk1_4+f0v1.*n2k1_2.*t6.*t8.*t98.*ve_v1.*xk1_2+f0d2.*t98.*t104.*t134.*t188.*ve_d2.*xk2_4),t187.*t209.*(-t104.*t134.*t188.*ve_d2.*xk2_4+f0v2.*n2k2_2.*t17.*t19.*t99.*ve_v2.*xk2_2+f0d3.*t99.*t107.*t132.*t190.*ve_d3.*xk3_4),t189.*t211.*(-t107.*t132.*t190.*ve_d3.*xk3_4+f0v3.*n2k3_2.*t25.*t27.*t100.*ve_v3.*xk3_2+f0d4.*t100.*t110.*t130.*t192.*ve_d4.*xk4_4),t191.*t213.*(-t110.*t130.*t192.*ve_d4.*xk4_4+f0v4.*n2k4_2.*t33.*t35.*t101.*ve_v4.*xk4_2+f0d5.*n2k5_2.*t43.*t101.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),t193.*t215.*(n2k5_2.*t41.*t43.*ve_v5.*xk5_2-n2k5_2.*t43.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0./xn3_4.^2.*(c33.*t11-c23.*xn3_3)-c33.*t4,t5.*t6.*(n2k1_3-n2k1_3.*t8.*ve_v1),t16.*t17.*(n2k2_3-n2k2_3.*t19.*ve_v2),t24.*t25.*(n2k3_3-n2k3_3.*t27.*ve_v3),t32.*t33.*(n2k4_3-n2k4_3.*t35.*ve_v4),t40.*t41.*(n2k5_3-n2k5_3.*t43.*ve_v5),t5.*t48.*(n2k1_3.*t49-n2k1_3.*t6.*t8.*ve_v1.*xk1_2),t16.*t58.*(n2k2_3.*t59-n2k2_3.*t17.*t19.*ve_v2.*xk2_2),t24.*t68.*(n2k3_3.*t69-n2k3_3.*t25.*t27.*ve_v3.*xk3_2),t32.*t78.*(n2k4_3.*t79-n2k4_3.*t33.*t35.*ve_v4.*xk4_2),t40.*t88.*(n2k5_3.*t89-n2k5_3.*t41.*t43.*ve_v5.*xk5_2),t121.*t122.*(t136+t146-t124.*t199.*ve_d1),t187.*t188.*(t137+t144-t104.*t145.*ve_d2),t189.*t190.*(t138+t142-t107.*t143.*ve_d3),t191.*t192.*(t139+t140-t110.*t141.*ve_d4),t193.*t194.*(f0v5.*n2k5_3.*t43.*t195.*ve_v5-n2k5_3.*t43.*t113.*ve_d5.*ve_v5),t121.*t197.*(-t122.*t124.*t199.*ve_d1.*xk1_4+f0v1.*n2k1_3.*t6.*t8.*t98.*ve_v1.*xk1_2+f0d2.*t98.*t104.*t145.*t188.*ve_d2.*xk2_4),t187.*t209.*(-t104.*t145.*t188.*ve_d2.*xk2_4+f0v2.*n2k2_3.*t17.*t19.*t99.*ve_v2.*xk2_2+f0d3.*t99.*t107.*t143.*t190.*ve_d3.*xk3_4),t189.*t211.*(-t107.*t143.*t190.*ve_d3.*xk3_4+f0v3.*n2k3_3.*t25.*t27.*t100.*ve_v3.*xk3_2+f0d4.*t100.*t110.*t141.*t192.*ve_d4.*xk4_4),t191.*t213.*(-t110.*t141.*t192.*ve_d4.*xk4_4+f0v4.*n2k4_3.*t33.*t35.*t101.*ve_v4.*xk4_2+f0d5.*n2k5_3.*t43.*t101.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),t193.*t215.*(n2k5_3.*t41.*t43.*ve_v5.*xk5_2-n2k5_3.*t43.*t113.*t194.*ve_d5.*ve_v5.*xk5_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5.*t50.*(t12+t13+t14-t8.*t55+1.0)-t6.*t8.*t15.*t57,0.0,0.0,0.0,0.0,t5.*t48.*(t8.*t50.*t55.*xk1_2-t6.*t8.*t15.*t57.*tt_v1.*xk1_2),0.0,0.0,0.0,0.0,t121.*t122.*(f0v1.*t8.*t15.*t57.*t98.*tt_v1-f0v1.*t8.*t15.*t57.*t98.*t124.*tt_v1.*ve_d1),0.0,0.0,0.0,0.0,-t121.*t197.*(f0v1.*t8.*t50.*t55.*t98.*xk1_2-f0v1.*t6.*t8.*t15.*t57.*t98.*tt_v1.*xk1_2+f0v1.*t8.*t15.*t57.*t98.*t122.*t124.*tt_v1.*ve_d1.*xk1_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t16.*t60.*(t20+t21+t22-t19.*t65+1.0)-t17.*t19.*t23.*t67,0.0,0.0,0.0,0.0,t16.*t58.*(t19.*t60.*t65.*xk2_2-t17.*t19.*t23.*t67.*tt_v2.*xk2_2),0.0,0.0,0.0,t121.*t122.*(f0v2.*t19.*t23.*t67.*t98.*t104.*tt_v2.*ve_d2-f0v2.*t19.*t23.*t67.*t98.*t104.*t124.*tt_v2.*ve_d1.*ve_d2),t187.*t188.*(f0v2.*t19.*t23.*t67.*t99.*tt_v2-f0v2.*t19.*t23.*t67.*t99.*t104.*tt_v2.*ve_d2),0.0,0.0,0.0,t121.*t197.*(f0v2.*t19.*t23.*t67.*t98.*t104.*t188.*tt_v2.*ve_d2.*xk2_4-f0v2.*t19.*t23.*t67.*t98.*t104.*t122.*t124.*tt_v2.*ve_d1.*ve_d2.*xk1_4),-t187.*t209.*(f0v2.*t19.*t60.*t65.*t99.*xk2_2-f0v2.*t17.*t19.*t23.*t67.*t99.*tt_v2.*xk2_2+f0v2.*t19.*t23.*t67.*t99.*t104.*t188.*tt_v2.*ve_d2.*xk2_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t24.*t70.*(t28+t29+t30-t27.*t75+1.0)-t25.*t27.*t31.*t77,0.0,0.0,0.0,0.0,t24.*t68.*(t27.*t70.*t75.*xk3_2-t25.*t27.*t31.*t77.*tt_v3.*xk3_2),0.0,0.0,t121.*t122.*(f0v3.*t27.*t31.*t77.*t98.*t104.*t107.*tt_v3.*ve_d2.*ve_d3-f0v3.*t27.*t31.*t77.*t98.*t104.*t107.*t124.*tt_v3.*ve_d1.*ve_d2.*ve_d3),t187.*t188.*(f0v3.*t27.*t31.*t77.*t99.*t107.*tt_v3.*ve_d3-f0v3.*t27.*t31.*t77.*t99.*t104.*t107.*tt_v3.*ve_d2.*ve_d3),t189.*t190.*(f0v3.*t27.*t31.*t77.*t100.*tt_v3-f0v3.*t27.*t31.*t77.*t100.*t107.*tt_v3.*ve_d3),0.0,0.0,t121.*t197.*(f0v3.*t27.*t31.*t77.*t98.*t104.*t107.*t188.*tt_v3.*ve_d2.*ve_d3.*xk2_4-f0v3.*t27.*t31.*t77.*t98.*t104.*t107.*t122.*t124.*tt_v3.*ve_d1.*ve_d2.*ve_d3.*xk1_4),t187.*t209.*(f0v3.*t27.*t31.*t77.*t99.*t107.*t190.*tt_v3.*ve_d3.*xk3_4-f0v3.*t27.*t31.*t77.*t99.*t104.*t107.*t188.*tt_v3.*ve_d2.*ve_d3.*xk2_4),-t189.*t211.*(f0v3.*t27.*t70.*t75.*t100.*xk3_2-f0v3.*t25.*t27.*t31.*t77.*t100.*tt_v3.*xk3_2+f0v3.*t27.*t31.*t77.*t100.*t107.*t190.*tt_v3.*ve_d3.*xk3_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t32.*t80.*(t36+t37+t38-t35.*t85+1.0)-t33.*t35.*t39.*t87,0.0,0.0,0.0,0.0,t32.*t78.*(t35.*t80.*t85.*xk4_2-t33.*t35.*t39.*t87.*tt_v4.*xk4_2),0.0,t121.*t122.*(f0v4.*t35.*t39.*t87.*t98.*t104.*t107.*t110.*tt_v4.*ve_d2.*ve_d3.*ve_d4-f0v4.*t35.*t39.*t87.*t98.*t104.*t107.*t110.*t124.*tt_v4.*ve_d1.*ve_d2.*ve_d3.*ve_d4),t187.*t188.*(f0v4.*t35.*t39.*t87.*t99.*t107.*t110.*tt_v4.*ve_d3.*ve_d4-f0v4.*t35.*t39.*t87.*t99.*t104.*t107.*t110.*tt_v4.*ve_d2.*ve_d3.*ve_d4),t189.*t190.*(f0v4.*t35.*t39.*t87.*t100.*t110.*tt_v4.*ve_d4-f0v4.*t35.*t39.*t87.*t100.*t107.*t110.*tt_v4.*ve_d3.*ve_d4),t191.*t192.*(f0v4.*t35.*t39.*t87.*t101.*tt_v4-f0v4.*t35.*t39.*t87.*t101.*t110.*tt_v4.*ve_d4),0.0,t121.*t197.*(f0v4.*t35.*t39.*t87.*t98.*t104.*t107.*t110.*t188.*tt_v4.*ve_d2.*ve_d3.*ve_d4.*xk2_4-f0v4.*t35.*t39.*t87.*t98.*t104.*t107.*t110.*t122.*t124.*tt_v4.*ve_d1.*ve_d2.*ve_d3.*ve_d4.*xk1_4),t187.*t209.*(f0v4.*t35.*t39.*t87.*t99.*t107.*t110.*t190.*tt_v4.*ve_d3.*ve_d4.*xk3_4-f0v4.*t35.*t39.*t87.*t99.*t104.*t107.*t110.*t188.*tt_v4.*ve_d2.*ve_d3.*ve_d4.*xk2_4),t189.*t211.*(f0v4.*t35.*t39.*t87.*t100.*t110.*t192.*tt_v4.*ve_d4.*xk4_4-f0v4.*t35.*t39.*t87.*t100.*t107.*t110.*t190.*tt_v4.*ve_d3.*ve_d4.*xk3_4),-t191.*t213.*(f0v4.*t35.*t80.*t85.*t101.*xk4_2-f0v4.*t33.*t35.*t39.*t87.*t101.*tt_v4.*xk4_2+f0v4.*t35.*t39.*t87.*t101.*t110.*t192.*tt_v4.*ve_d4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t40.*t90.*(t44+t45+t46-t43.*t95+1.0)-t41.*t43.*t47.*t97,0.0,0.0,0.0,0.0,t40.*t88.*(t216-t41.*t43.*t47.*t97.*tt_v5.*xk5_2),t121.*t122.*(f0d5.*t43.*t47.*t97.*t98.*t104.*t107.*t110.*t113.*tt_v5.*ve_d2.*ve_d3.*ve_d4.*ve_d5-f0d5.*t43.*t47.*t97.*t98.*t104.*t107.*t110.*t113.*t124.*tt_v5.*ve_d1.*ve_d2.*ve_d3.*ve_d4.*ve_d5),t187.*t188.*(f0d5.*t43.*t47.*t97.*t99.*t107.*t110.*t113.*tt_v5.*ve_d3.*ve_d4.*ve_d5-f0d5.*t43.*t47.*t97.*t99.*t104.*t107.*t110.*t113.*tt_v5.*ve_d2.*ve_d3.*ve_d4.*ve_d5),t189.*t190.*(f0d5.*t43.*t47.*t97.*t100.*t110.*t113.*tt_v5.*ve_d4.*ve_d5-f0d5.*t43.*t47.*t97.*t100.*t107.*t110.*t113.*tt_v5.*ve_d3.*ve_d4.*ve_d5),t191.*t192.*(f0d5.*t43.*t47.*t97.*t101.*t113.*tt_v5.*ve_d5-f0d5.*t43.*t47.*t97.*t101.*t110.*t113.*tt_v5.*ve_d4.*ve_d5),t193.*t194.*(f0v5.*t43.*t47.*t97.*t195.*tt_v5-t43.*t47.*t97.*t113.*tt_v5.*ve_d5),t121.*t197.*(f0d5.*t43.*t47.*t97.*t98.*t104.*t107.*t110.*t113.*t188.*tt_v5.*ve_d2.*ve_d3.*ve_d4.*ve_d5.*xk2_4-f0d5.*t43.*t47.*t97.*t98.*t104.*t107.*t110.*t113.*t122.*t124.*tt_v5.*ve_d1.*ve_d2.*ve_d3.*ve_d4.*ve_d5.*xk1_4),t187.*t209.*(f0d5.*t43.*t47.*t97.*t99.*t107.*t110.*t113.*t190.*tt_v5.*ve_d3.*ve_d4.*ve_d5.*xk3_4-f0d5.*t43.*t47.*t97.*t99.*t104.*t107.*t110.*t113.*t188.*tt_v5.*ve_d2.*ve_d3.*ve_d4.*ve_d5.*xk2_4),t189.*t211.*(f0d5.*t43.*t47.*t97.*t100.*t110.*t113.*t192.*tt_v5.*ve_d4.*ve_d5.*xk4_4-f0d5.*t43.*t47.*t97.*t100.*t107.*t110.*t113.*t190.*tt_v5.*ve_d3.*ve_d4.*ve_d5.*xk3_4),t191.*t213.*(f0d5.*t43.*t47.*t97.*t101.*t113.*t194.*tt_v5.*ve_d5.*xk5_4-f0d5.*t43.*t47.*t97.*t101.*t110.*t113.*t192.*tt_v5.*ve_d4.*ve_d5.*xk4_4),-t193.*t215.*(t216-t41.*t43.*t47.*t97.*tt_v5.*xk5_2+t43.*t47.*t97.*t113.*t194.*tt_v5.*ve_d5.*xk5_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5.*1.0./xk1_2.^2.*(t49.*(nr1+t12+t13+t14)-t6.*t8.*t55.*xk1_2)-t5.*t6.*t8.*t48.*t55,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v1.*t6.*t8.*t55.*t98.*t121.*t197,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t16.*1.0./xk2_2.^2.*(t59.*(nr2+t20+t21+t22)-t17.*t19.*t65.*xk2_2)-t16.*t17.*t19.*t58.*t65,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v2.*t17.*t19.*t65.*t99.*t187.*t209,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t24.*1.0./xk3_2.^2.*(t69.*(nr3+t28+t29+t30)-t25.*t27.*t75.*xk3_2)-t24.*t25.*t27.*t68.*t75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v3.*t25.*t27.*t75.*t100.*t189.*t211,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t32.*1.0./xk4_2.^2.*(t79.*(nr4+t36+t37+t38)-t33.*t35.*t85.*xk4_2)-t32.*t33.*t35.*t78.*t85,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0v4.*t33.*t35.*t85.*t101.*t191.*t213,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t40.*1.0./xk5_2.^2.*(t89.*(nr5+t44+t45+t46)-t41.*t43.*t95.*xk5_2)-t40.*t41.*t43.*t88.*t95,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t41.*t43.*t95.*t193.*t215,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t121.*t200.*(t147+t177-t124.*t205)-t122.*t124.*t178.*t207,0.0,0.0,0.0,0.0,t121.*t197.*(t124.*t200.*t205.*xk1_4-t122.*t124.*t178.*t207.*tt_d1.*xk1_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t121.*t122.*(f0d2.*t98.*t104.*t148.*t180.*tt_d2-f0d2.*t98.*t104.*t124.*t148.*t180.*tt_d2.*ve_d1),-t187.*t208.*(t151+t173-t104.*t176)-t104.*t148.*t180.*t188,0.0,0.0,0.0,-t121.*t197.*(f0d2.*t98.*t104.*t176.*t208.*xk2_4-f0d2.*t98.*t104.*t148.*t180.*t188.*tt_d2.*xk2_4+f0d2.*t98.*t104.*t122.*t124.*t148.*t180.*tt_d2.*ve_d1.*xk1_4),t187.*t209.*(t104.*t176.*t208.*xk2_4-t104.*t148.*t180.*t188.*tt_d2.*xk2_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t121.*t122.*(f0d3.*t98.*t104.*t107.*t169.*t182.*tt_d3.*ve_d2-f0d3.*t98.*t104.*t107.*t124.*t169.*t182.*tt_d3.*ve_d1.*ve_d2),t187.*t188.*(f0d3.*t99.*t107.*t169.*t182.*tt_d3-f0d3.*t99.*t104.*t107.*t169.*t182.*tt_d3.*ve_d2),-t189.*t210.*(t152+t166-t107.*t172)-t107.*t169.*t182.*t190,0.0,0.0,t121.*t197.*(f0d3.*t98.*t104.*t107.*t169.*t182.*t188.*tt_d3.*ve_d2.*xk2_4-f0d3.*t98.*t104.*t107.*t122.*t124.*t169.*t182.*tt_d3.*ve_d1.*ve_d2.*xk1_4),-t187.*t209.*(f0d3.*t99.*t107.*t172.*t210.*xk3_4-f0d3.*t99.*t107.*t169.*t182.*t190.*tt_d3.*xk3_4+f0d3.*t99.*t104.*t107.*t169.*t182.*t188.*tt_d3.*ve_d2.*xk2_4),t189.*t211.*(t107.*t172.*t210.*xk3_4-t107.*t169.*t182.*t190.*tt_d3.*xk3_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t121.*t122.*(f0d4.*t98.*t104.*t107.*t110.*t153.*t184.*tt_d4.*ve_d2.*ve_d3-f0d4.*t98.*t104.*t107.*t110.*t124.*t153.*t184.*tt_d4.*ve_d1.*ve_d2.*ve_d3),t187.*t188.*(f0d4.*t99.*t107.*t110.*t153.*t184.*tt_d4.*ve_d3-f0d4.*t99.*t104.*t107.*t110.*t153.*t184.*tt_d4.*ve_d2.*ve_d3),t189.*t190.*(f0d4.*t100.*t110.*t153.*t184.*tt_d4-f0d4.*t100.*t107.*t110.*t153.*t184.*tt_d4.*ve_d3),-t191.*t212.*(t156+t162-t110.*t165)-t110.*t153.*t184.*t192,0.0,t121.*t197.*(f0d4.*t98.*t104.*t107.*t110.*t153.*t184.*t188.*tt_d4.*ve_d2.*ve_d3.*xk2_4-f0d4.*t98.*t104.*t107.*t110.*t122.*t124.*t153.*t184.*tt_d4.*ve_d1.*ve_d2.*ve_d3.*xk1_4),t187.*t209.*(f0d4.*t99.*t107.*t110.*t153.*t184.*t190.*tt_d4.*ve_d3.*xk3_4-f0d4.*t99.*t104.*t107.*t110.*t153.*t184.*t188.*tt_d4.*ve_d2.*ve_d3.*xk2_4),-t189.*t211.*(f0d4.*t100.*t110.*t165.*t212.*xk4_4-f0d4.*t100.*t110.*t153.*t184.*t192.*tt_d4.*xk4_4+f0d4.*t100.*t107.*t110.*t153.*t184.*t190.*tt_d4.*ve_d3.*xk3_4),t191.*t213.*(t110.*t165.*t212.*xk4_4-t110.*t153.*t184.*t192.*tt_d4.*xk4_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t121.*t122.*(f0d5.*t98.*t104.*t107.*t110.*t113.*t157.*t186.*tt_d5.*ve_d2.*ve_d3.*ve_d4-f0d5.*t98.*t104.*t107.*t110.*t113.*t124.*t157.*t186.*tt_d5.*ve_d1.*ve_d2.*ve_d3.*ve_d4),t187.*t188.*(f0d5.*t99.*t107.*t110.*t113.*t157.*t186.*tt_d5.*ve_d3.*ve_d4-f0d5.*t99.*t104.*t107.*t110.*t113.*t157.*t186.*tt_d5.*ve_d2.*ve_d3.*ve_d4),t189.*t190.*(f0d5.*t100.*t110.*t113.*t157.*t186.*tt_d5.*ve_d4-f0d5.*t100.*t107.*t110.*t113.*t157.*t186.*tt_d5.*ve_d3.*ve_d4),t191.*t192.*(f0d5.*t101.*t113.*t157.*t186.*tt_d5-f0d5.*t101.*t110.*t113.*t157.*t186.*tt_d5.*ve_d4),t193.*t214.*(t113.*t161-f0v5.*t43.*t95.*t195)-t113.*t157.*t186.*t194,t121.*t197.*(f0d5.*t98.*t104.*t107.*t110.*t113.*t157.*t186.*t188.*tt_d5.*ve_d2.*ve_d3.*ve_d4.*xk2_4-f0d5.*t98.*t104.*t107.*t110.*t113.*t122.*t124.*t157.*t186.*tt_d5.*ve_d1.*ve_d2.*ve_d3.*ve_d4.*xk1_4),t187.*t209.*(f0d5.*t99.*t107.*t110.*t113.*t157.*t186.*t190.*tt_d5.*ve_d3.*ve_d4.*xk3_4-f0d5.*t99.*t104.*t107.*t110.*t113.*t157.*t186.*t188.*tt_d5.*ve_d2.*ve_d3.*ve_d4.*xk2_4),t189.*t211.*(f0d5.*t100.*t110.*t113.*t157.*t186.*t192.*tt_d5.*ve_d4.*xk4_4-f0d5.*t100.*t107.*t110.*t113.*t157.*t186.*t190.*tt_d5.*ve_d3.*ve_d4.*xk3_4),-t191.*t213.*(f0d5.*t101.*t113.*t161.*t214.*xk5_4-f0d5.*t101.*t113.*t157.*t186.*t194.*tt_d5.*xk5_4+f0d5.*t101.*t110.*t113.*t157.*t186.*t192.*tt_d5.*ve_d4.*xk4_4),t193.*t215.*(t113.*t161.*t214.*xk5_4-t113.*t157.*t186.*t194.*tt_d5.*xk5_4),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t121.*1.0./xk1_4.^2.*(-t122.*t124.*t205.*xk1_4+f0d2.*t98.*t104.*t176.*t188.*xk2_4+f0v1.*t6.*t8.*t55.*t98.*xk1_2)-t121.*t122.*t124.*t197.*t205,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d2.*t98.*t104.*t121.*t176.*t188.*t197,-t187.*1.0./xk2_4.^2.*(-t104.*t176.*t188.*xk2_4+f0d3.*t99.*t107.*t172.*t190.*xk3_4+f0v2.*t17.*t19.*t65.*t99.*xk2_2)-t104.*t176.*t187.*t188.*t209,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d3.*t99.*t107.*t172.*t187.*t190.*t209,-t189.*1.0./xk3_4.^2.*(-t107.*t172.*t190.*xk3_4+f0d4.*t100.*t110.*t165.*t192.*xk4_4+f0v3.*t25.*t27.*t75.*t100.*xk3_2)-t107.*t172.*t189.*t190.*t211,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d4.*t100.*t110.*t165.*t189.*t192.*t211,-t191.*1.0./xk4_4.^2.*(-t110.*t165.*t192.*xk4_4+f0d5.*t101.*t113.*t161.*t194.*xk5_4+f0v4.*t33.*t35.*t85.*t101.*xk4_2)-t110.*t165.*t191.*t192.*t213,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f0d5.*t101.*t113.*t161.*t191.*t194.*t213,-t193.*1.0./xk5_4.^2.*(t41.*t43.*t95.*xk5_2-t113.*t161.*t194.*xk5_4)-t113.*t161.*t193.*t194.*t215],[32,32]);
