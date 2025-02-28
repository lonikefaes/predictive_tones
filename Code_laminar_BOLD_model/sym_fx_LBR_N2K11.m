function f_s = sym_fx_LBR_N2K11(xn1_1,xn2_1,xn1_2,xn2_2,xn1_3,xn2_3,xn1_4,xn2_4,xk1_1,xk2_1,xk3_1,xk4_1,xk5_1,xk6_1,xk7_1,xk8_1,xk9_1,xk10_1,xk11_1,xk1_2,xk2_2,xk3_2,xk4_2,xk5_2,xk6_2,xk7_2,xk8_2,xk9_2,xk10_2,xk11_2,xk1_3,xk2_3,xk3_3,xk4_3,xk5_3,xk6_3,xk7_3,xk8_3,xk9_3,xk10_3,xk11_3,xk1_4,xk2_4,xk3_4,xk4_4,xk5_4,xk6_4,xk7_4,xk8_4,xk9_4,xk10_4,xk11_4,a1_1,a2_1,a1_2,a2_2,mu1,mu2,cu1,cu2,lam1,lam2,c11,c12,c21,c22,c31,c32,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k5_1,n2k6_1,n2k7_1,n2k8_1,n2k9_1,n2k10_1,n2k11_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k5_2,n2k6_2,n2k7_2,n2k8_2,n2k9_2,n2k10_2,n2k11_2,al_v1,al_v2,al_v3,al_v4,al_v5,al_v6,al_v7,al_v8,al_v9,al_v10,al_v11,al_d1,al_d2,al_d3,al_d4,al_d5,al_d6,al_d7,al_d8,al_d9,al_d10,al_d11,tt_v1,tt_v2,tt_v3,tt_v4,tt_v5,tt_v6,tt_v7,tt_v8,tt_v9,tt_v10,tt_v11,tt_d1,tt_d2,tt_d3,tt_d4,tt_d5,tt_d6,tt_d7,tt_d8,tt_d9,tt_d10,tt_d11,nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,ve_v1,ve_v2,ve_v3,ve_v4,ve_v5,ve_v6,ve_v7,ve_v8,ve_v9,ve_v10,ve_v11,ve_d1,ve_d2,ve_d3,ve_d4,ve_d5,ve_d6,ve_d7,ve_d8,ve_d9,ve_d10,ve_d11,f0v1,f0v2,f0v3,f0v4,f0v5,f0v6,f0v7,f0v8,f0v9,f0v10,f0v11,f0d1,f0d2,f0d3,f0d4,f0d5,f0d6,f0d7,f0d8,f0d9,f0d10,f0d11)
%SYM_FX_LBR_N2K11
%    F_S = SYM_FX_LBR_N2K11(XN1_1,XN2_1,XN1_2,XN2_2,XN1_3,XN2_3,XN1_4,XN2_4,XK1_1,XK2_1,XK3_1,XK4_1,XK5_1,XK6_1,XK7_1,XK8_1,XK9_1,XK10_1,XK11_1,XK1_2,XK2_2,XK3_2,XK4_2,XK5_2,XK6_2,XK7_2,XK8_2,XK9_2,XK10_2,XK11_2,XK1_3,XK2_3,XK3_3,XK4_3,XK5_3,XK6_3,XK7_3,XK8_3,XK9_3,XK10_3,XK11_3,XK1_4,XK2_4,XK3_4,XK4_4,XK5_4,XK6_4,XK7_4,XK8_4,XK9_4,XK10_4,XK11_4,A1_1,A2_1,A1_2,A2_2,MU1,MU2,CU1,CU2,LAM1,LAM2,C11,C12,C21,C22,C31,C32,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K5_1,N2K6_1,N2K7_1,N2K8_1,N2K9_1,N2K10_1,N2K11_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K5_2,N2K6_2,N2K7_2,N2K8_2,N2K9_2,N2K10_2,N2K11_2,AL_V1,AL_V2,AL_V3,AL_V4,AL_V5,AL_V6,AL_V7,AL_V8,AL_V9,AL_V10,AL_V11,AL_D1,AL_D2,AL_D3,AL_D4,AL_D5,AL_D6,AL_D7,AL_D8,AL_D9,AL_D10,AL_D11,TT_V1,TT_V2,TT_V3,TT_V4,TT_V5,TT_V6,TT_V7,TT_V8,TT_V9,TT_V10,TT_V11,TT_D1,TT_D2,TT_D3,TT_D4,TT_D5,TT_D6,TT_D7,TT_D8,TT_D9,TT_D10,TT_D11,NR1,NR2,NR3,NR4,NR5,NR6,NR7,NR8,NR9,NR10,NR11,VE_V1,VE_V2,VE_V3,VE_V4,VE_V5,VE_V6,VE_V7,VE_V8,VE_V9,VE_V10,VE_V11,VE_D1,VE_D2,VE_D3,VE_D4,VE_D5,VE_D6,VE_D7,VE_D8,VE_D9,VE_D10,VE_D11,F0V1,F0V2,F0V3,F0V4,F0V5,F0V6,F0V7,F0V8,F0V9,F0V10,F0V11,F0D1,F0D2,F0D3,F0D4,F0D5,F0D6,F0D7,F0D8,F0D9,F0D10,F0D11)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    06-Feb-2020 13:35:55

t2 = xn1_4-1.0;
t3 = xn2_4-1.0;
t4 = n2k1_1.*t2;
t5 = n2k1_2.*t3;
t6 = n2k2_1.*t2;
t7 = n2k2_2.*t3;
t8 = n2k3_1.*t2;
t9 = n2k3_2.*t3;
t10 = n2k4_1.*t2;
t11 = n2k4_2.*t3;
t12 = n2k5_1.*t2;
t13 = n2k5_2.*t3;
t14 = n2k6_1.*t2;
t15 = n2k6_2.*t3;
t16 = n2k7_1.*t2;
t17 = n2k7_2.*t3;
t18 = n2k8_1.*t2;
t19 = n2k8_2.*t3;
t20 = n2k9_1.*t2;
t21 = n2k9_2.*t3;
t22 = n2k10_1.*t2;
t23 = n2k10_2.*t3;
t24 = n2k11_1.*t2;
t25 = n2k11_2.*t3;
t26 = 1.0./tt_v1;
t27 = 1.0./xk1_1;
t28 = tt_v1+ve_v1;
t29 = 1.0./t28;
t30 = 1.0./al_v1;
t31 = xk1_1.^t30;
t32 = t31.*tt_v1;
t33 = t4+t5+1.0;
t34 = t33.*ve_v1;
t35 = t32+t34;
t36 = 1.0./tt_v2;
t37 = 1.0./xk2_1;
t38 = tt_v2+ve_v2;
t39 = 1.0./t38;
t40 = 1.0./al_v2;
t41 = xk2_1.^t40;
t42 = t41.*tt_v2;
t43 = t6+t7+1.0;
t44 = t43.*ve_v2;
t45 = t42+t44;
t46 = 1.0./tt_v3;
t47 = 1.0./xk3_1;
t48 = tt_v3+ve_v3;
t49 = 1.0./t48;
t50 = 1.0./al_v3;
t51 = xk3_1.^t50;
t52 = t51.*tt_v3;
t53 = t8+t9+1.0;
t54 = t53.*ve_v3;
t55 = t52+t54;
t56 = 1.0./tt_v4;
t57 = 1.0./xk4_1;
t58 = tt_v4+ve_v4;
t59 = 1.0./t58;
t60 = 1.0./al_v4;
t61 = xk4_1.^t60;
t62 = t61.*tt_v4;
t63 = t10+t11+1.0;
t64 = t63.*ve_v4;
t65 = t62+t64;
t66 = 1.0./tt_v5;
t67 = 1.0./xk5_1;
t68 = tt_v5+ve_v5;
t69 = 1.0./t68;
t70 = 1.0./al_v5;
t71 = xk5_1.^t70;
t72 = t71.*tt_v5;
t73 = t12+t13+1.0;
t74 = t73.*ve_v5;
t75 = t72+t74;
t76 = 1.0./tt_v6;
t77 = 1.0./xk6_1;
t78 = tt_v6+ve_v6;
t79 = 1.0./t78;
t80 = 1.0./al_v6;
t81 = xk6_1.^t80;
t82 = t81.*tt_v6;
t83 = t14+t15+1.0;
t84 = t83.*ve_v6;
t85 = t82+t84;
t86 = 1.0./tt_v7;
t87 = 1.0./xk7_1;
t88 = tt_v7+ve_v7;
t89 = 1.0./t88;
t90 = 1.0./al_v7;
t91 = xk7_1.^t90;
t92 = t91.*tt_v7;
t93 = t16+t17+1.0;
t94 = t93.*ve_v7;
t95 = t92+t94;
t96 = 1.0./tt_v8;
t97 = 1.0./xk8_1;
t98 = tt_v8+ve_v8;
t99 = 1.0./t98;
t100 = 1.0./al_v8;
t101 = xk8_1.^t100;
t102 = t101.*tt_v8;
t103 = t18+t19+1.0;
t104 = t103.*ve_v8;
t105 = t102+t104;
t106 = 1.0./tt_v9;
t107 = 1.0./xk9_1;
t108 = tt_v9+ve_v9;
t109 = 1.0./t108;
t110 = 1.0./al_v9;
t111 = xk9_1.^t110;
t112 = t111.*tt_v9;
t113 = t20+t21+1.0;
t114 = t113.*ve_v9;
t115 = t112+t114;
t116 = 1.0./tt_v10;
t117 = 1.0./xk10_1;
t118 = tt_v10+ve_v10;
t119 = 1.0./t118;
t120 = 1.0./al_v10;
t121 = xk10_1.^t120;
t122 = t121.*tt_v10;
t123 = t22+t23+1.0;
t124 = t123.*ve_v10;
t125 = t122+t124;
t126 = 1.0./tt_v11;
t127 = 1.0./xk11_1;
t128 = tt_v11+ve_v11;
t129 = 1.0./t128;
t130 = 1.0./al_v11;
t131 = xk11_1.^t130;
t132 = t131.*tt_v11;
t133 = t24+t25+1.0;
t134 = t133.*ve_v11;
t135 = t132+t134;
t136 = 1.0./f0d3;
t137 = 1.0./f0d6;
t138 = 1.0./f0d7;
t139 = 1.0./f0d9;
t140 = 1.0./f0d10;
t141 = 1.0./f0d8;
t142 = 1.0./f0d5;
t143 = 1.0./f0d4;
t144 = 1.0./f0d2;
t145 = 1.0./f0d1;
t146 = tt_d3+ve_d3;
t147 = 1.0./t146;
t148 = 1.0./al_d3;
t149 = xk3_3.^t148;
t150 = t149.*tt_d3;
t151 = f0v3.*t49.*t55.*t136;
t152 = tt_d4+ve_d4;
t153 = 1.0./t152;
t154 = 1.0./al_d4;
t155 = xk4_3.^t154;
t156 = t155.*tt_d4;
t157 = tt_d5+ve_d5;
t158 = 1.0./t157;
t159 = f0v6.*t79.*t85.*t137;
t160 = tt_d7+ve_d7;
t161 = 1.0./t160;
t162 = 1.0./al_d7;
t163 = xk7_3.^t162;
t164 = t163.*tt_d7;
t165 = f0v7.*t89.*t95.*t138;
t166 = tt_d8+ve_d8;
t167 = 1.0./t166;
t168 = 1.0./al_d8;
t169 = xk8_3.^t168;
t170 = t169.*tt_d8;
t171 = f0v9.*t109.*t115.*t139;
t172 = tt_d10+ve_d10;
t173 = 1.0./t172;
t174 = 1.0./al_d10;
t175 = xk10_3.^t174;
t176 = t175.*tt_d10;
t177 = f0v10.*t119.*t125.*t140;
t178 = tt_d11+ve_d11;
t179 = 1.0./t178;
t180 = 1.0./al_d11;
t181 = xk11_3.^t180;
t182 = t181.*tt_d11;
t183 = t129.*t135.*ve_d11;
t184 = t182+t183;
t185 = f0d11.*t140.*t179.*t184;
t186 = t177+t185;
t187 = t186.*ve_d10;
t188 = t176+t187;
t189 = f0d10.*t139.*t173.*t188;
t190 = t171+t189;
t191 = t190.*ve_d9;
t192 = 1.0./al_d9;
t193 = xk9_3.^t192;
t194 = t193.*tt_d9;
t195 = t191+t194;
t196 = tt_d9+ve_d9;
t197 = 1.0./t196;
t198 = f0d9.*t141.*t195.*t197;
t199 = f0v8.*t99.*t105.*t141;
t200 = t198+t199;
t201 = t200.*ve_d8;
t202 = t170+t201;
t203 = f0d8.*t138.*t167.*t202;
t204 = t165+t203;
t205 = t204.*ve_d7;
t206 = t164+t205;
t207 = f0d7.*t137.*t161.*t206;
t208 = t159+t207;
t209 = t208.*ve_d6;
t210 = 1.0./al_d6;
t211 = xk6_3.^t210;
t212 = t211.*tt_d6;
t213 = t209+t212;
t214 = tt_d6+ve_d6;
t215 = 1.0./t214;
t216 = f0d6.*t142.*t213.*t215;
t217 = f0v5.*t69.*t75.*t142;
t218 = t216+t217;
t219 = t218.*ve_d5;
t220 = 1.0./al_d5;
t221 = xk5_3.^t220;
t222 = t221.*tt_d5;
t223 = t219+t222;
t224 = f0d5.*t143.*t158.*t223;
t225 = f0v4.*t59.*t65.*t143;
t226 = t224+t225;
t227 = t226.*ve_d4;
t228 = t156+t227;
t229 = f0d4.*t136.*t153.*t228;
t230 = t151+t229;
t231 = t230.*ve_d3;
t232 = t150+t231;
t233 = f0d3.*t144.*t147.*t232;
t234 = f0v2.*t39.*t45.*t144;
t235 = t233+t234;
t236 = t235.*ve_d2;
t237 = 1.0./al_d2;
t238 = xk2_3.^t237;
t239 = t238.*tt_d2;
t240 = t236+t239;
t241 = tt_d2+ve_d2;
t242 = 1.0./t241;
t243 = f0d2.*t145.*t240.*t242;
t244 = f0v1.*t29.*t35.*t145;
t245 = 1.0./tt_d1;
t246 = tt_d1+ve_d1;
t247 = 1.0./t246;
t248 = 1.0./al_d1;
t249 = xk1_3.^t248;
t250 = t249.*tt_d1;
t251 = t243+t244;
t252 = t251.*ve_d1;
t253 = t250+t252;
t254 = 1.0./tt_d2;
t255 = 1.0./xk2_3;
t256 = 1.0./tt_d3;
t257 = 1.0./xk3_3;
t258 = 1.0./tt_d4;
t259 = 1.0./xk4_3;
t260 = 1.0./tt_d5;
t261 = 1.0./xk5_3;
t262 = 1.0./tt_d6;
t263 = 1.0./xk6_3;
t264 = 1.0./tt_d7;
t265 = 1.0./xk7_3;
t266 = 1.0./tt_d8;
t267 = 1.0./xk8_3;
t268 = 1.0./tt_d9;
t269 = 1.0./xk9_3;
t270 = 1.0./tt_d10;
t271 = 1.0./xk10_3;
t272 = 1.0./xk11_3;
t273 = 1.0./tt_d11;
f_s = [cu1+a1_1.*xn1_1+a1_2.*xn2_1-mu1.*xn1_2;cu2+a2_1.*xn1_1+a2_2.*xn2_1-mu2.*xn2_2;lam1.*(xn1_1-xn1_2);lam2.*(xn2_1-xn2_2);xn1_1-c11.*xn1_3;xn2_1-c12.*xn2_3;-(c31.*t2-c21.*xn1_3)./xn1_4;-(c32.*t3-c22.*xn2_3)./xn2_4;t26.*t27.*(t4+t5-t29.*t35+1.0);t36.*t37.*(t6+t7-t39.*t45+1.0);t46.*t47.*(t8+t9-t49.*t55+1.0);t56.*t57.*(t10+t11-t59.*t65+1.0);t66.*t67.*(t12+t13-t69.*t75+1.0);t76.*t77.*(t14+t15-t79.*t85+1.0);t86.*t87.*(t16+t17-t89.*t95+1.0);t96.*t97.*(t18+t19-t99.*t105+1.0);t106.*t107.*(t20+t21-t109.*t115+1.0);t116.*t117.*(t22+t23-t119.*t125+1.0);t126.*t127.*(t24+t25-t129.*t135+1.0);(t26.*((nr1+t4+t5)./nr1-t27.*t29.*t35.*xk1_2))./xk1_2;(t36.*((nr2+t6+t7)./nr2-t37.*t39.*t45.*xk2_2))./xk2_2;(t46.*((nr3+t8+t9)./nr3-t47.*t49.*t55.*xk3_2))./xk3_2;(t56.*((nr4+t10+t11)./nr4-t57.*t59.*t65.*xk4_2))./xk4_2;(t66.*((nr5+t12+t13)./nr5-t67.*t69.*t75.*xk5_2))./xk5_2;(t76.*((nr6+t14+t15)./nr6-t77.*t79.*t85.*xk6_2))./xk6_2;(t86.*((nr7+t16+t17)./nr7-t87.*t89.*t95.*xk7_2))./xk7_2;(t96.*((nr8+t18+t19)./nr8-t97.*t99.*t105.*xk8_2))./xk8_2;(t106.*((nr9+t20+t21)./nr9-t107.*t109.*t115.*xk9_2))./xk9_2;(t116.*((nr10+t22+t23)./nr10-t117.*t119.*t125.*xk10_2))./xk10_2;(t126.*((nr11+t24+t25)./nr11-t127.*t129.*t135.*xk11_2))./xk11_2;(t245.*(t243+t244-t247.*t253))./ve_d1;(t254.*(t233+t234-t240.*t242))./ve_d2;(t256.*(t151+t229-t147.*t232))./ve_d3;(t258.*(t224+t225-t153.*t228))./ve_d4;(t260.*(t216+t217-t158.*t223))./ve_d5;(t262.*(t159+t207-t213.*t215))./ve_d6;(t264.*(t165+t203-t161.*t206))./ve_d7;(t266.*(t198+t199-t167.*t202))./ve_d8;(t268.*(t171+t189-t195.*t197))./ve_d9;(t270.*(t177+t185-t173.*t188))./ve_d10;-t272.*t273.*(t179.*t184-(f0v11.*t129.*t135)./f0d11);(t245.*(-(t247.*t253.*xk1_4)./xk1_3+f0d2.*t145.*t240.*t242.*t255.*xk2_4+f0v1.*t27.*t29.*t35.*t145.*xk1_2))./xk1_4;(t254.*(-t240.*t242.*t255.*xk2_4+f0d3.*t144.*t147.*t232.*t257.*xk3_4+f0v2.*t37.*t39.*t45.*t144.*xk2_2))./xk2_4;(t256.*(-t147.*t232.*t257.*xk3_4+f0d4.*t136.*t153.*t228.*t259.*xk4_4+f0v3.*t47.*t49.*t55.*t136.*xk3_2))./xk3_4;(t258.*(-t153.*t228.*t259.*xk4_4+f0d5.*t143.*t158.*t223.*t261.*xk5_4+f0v4.*t57.*t59.*t65.*t143.*xk4_2))./xk4_4;(t260.*(-t158.*t223.*t261.*xk5_4+f0d6.*t142.*t213.*t215.*t263.*xk6_4+f0v5.*t67.*t69.*t75.*t142.*xk5_2))./xk5_4;(t262.*(-t213.*t215.*t263.*xk6_4+f0d7.*t137.*t161.*t206.*t265.*xk7_4+f0v6.*t77.*t79.*t85.*t137.*xk6_2))./xk6_4;(t264.*(-t161.*t206.*t265.*xk7_4+f0d8.*t138.*t167.*t202.*t267.*xk8_4+f0v7.*t87.*t89.*t95.*t138.*xk7_2))./xk7_4;(t266.*(-t167.*t202.*t267.*xk8_4+f0d9.*t141.*t195.*t197.*t269.*xk9_4+f0v8.*t97.*t99.*t105.*t141.*xk8_2))./xk8_4;(t268.*(-t195.*t197.*t269.*xk9_4+f0d10.*t139.*t173.*t188.*t271.*xk10_4+f0v9.*t107.*t109.*t115.*t139.*xk9_2))./xk9_4;(t270.*(-t173.*t188.*t271.*xk10_4+f0d11.*t140.*t179.*t184.*t272.*xk11_4+f0v10.*t117.*t119.*t125.*t140.*xk10_2))./xk10_4;(t273.*(t127.*t129.*t135.*xk11_2-t179.*t184.*t272.*xk11_4))./xk11_4];
