function f_s = sym_fx_LBR_N2K12(xn1_1,xn2_1,xn1_2,xn2_2,xn1_3,xn2_3,xn1_4,xn2_4,xk1_1,xk2_1,xk3_1,xk4_1,xk5_1,xk6_1,xk7_1,xk8_1,xk9_1,xk10_1,xk11_1,xk12_1,xk1_2,xk2_2,xk3_2,xk4_2,xk5_2,xk6_2,xk7_2,xk8_2,xk9_2,xk10_2,xk11_2,xk12_2,xk1_3,xk2_3,xk3_3,xk4_3,xk5_3,xk6_3,xk7_3,xk8_3,xk9_3,xk10_3,xk11_3,xk12_3,xk1_4,xk2_4,xk3_4,xk4_4,xk5_4,xk6_4,xk7_4,xk8_4,xk9_4,xk10_4,xk11_4,xk12_4,a1_1,a2_1,a1_2,a2_2,mu1,mu2,cu1,cu2,lam1,lam2,c11,c12,c21,c22,c31,c32,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k5_1,n2k6_1,n2k7_1,n2k8_1,n2k9_1,n2k10_1,n2k11_1,n2k12_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k5_2,n2k6_2,n2k7_2,n2k8_2,n2k9_2,n2k10_2,n2k11_2,n2k12_2,al_v1,al_v2,al_v3,al_v4,al_v5,al_v6,al_v7,al_v8,al_v9,al_v10,al_v11,al_v12,al_d1,al_d2,al_d3,al_d4,al_d5,al_d6,al_d7,al_d8,al_d9,al_d10,al_d11,al_d12,tt_v1,tt_v2,tt_v3,tt_v4,tt_v5,tt_v6,tt_v7,tt_v8,tt_v9,tt_v10,tt_v11,tt_v12,tt_d1,tt_d2,tt_d3,tt_d4,tt_d5,tt_d6,tt_d7,tt_d8,tt_d9,tt_d10,tt_d11,tt_d12,nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,ve_v1,ve_v2,ve_v3,ve_v4,ve_v5,ve_v6,ve_v7,ve_v8,ve_v9,ve_v10,ve_v11,ve_v12,ve_d1,ve_d2,ve_d3,ve_d4,ve_d5,ve_d6,ve_d7,ve_d8,ve_d9,ve_d10,ve_d11,ve_d12,f0v1,f0v2,f0v3,f0v4,f0v5,f0v6,f0v7,f0v8,f0v9,f0v10,f0v11,f0v12,f0d1,f0d2,f0d3,f0d4,f0d5,f0d6,f0d7,f0d8,f0d9,f0d10,f0d11,f0d12)
%SYM_FX_LBR_N2K12
%    F_S = SYM_FX_LBR_N2K12(XN1_1,XN2_1,XN1_2,XN2_2,XN1_3,XN2_3,XN1_4,XN2_4,XK1_1,XK2_1,XK3_1,XK4_1,XK5_1,XK6_1,XK7_1,XK8_1,XK9_1,XK10_1,XK11_1,XK12_1,XK1_2,XK2_2,XK3_2,XK4_2,XK5_2,XK6_2,XK7_2,XK8_2,XK9_2,XK10_2,XK11_2,XK12_2,XK1_3,XK2_3,XK3_3,XK4_3,XK5_3,XK6_3,XK7_3,XK8_3,XK9_3,XK10_3,XK11_3,XK12_3,XK1_4,XK2_4,XK3_4,XK4_4,XK5_4,XK6_4,XK7_4,XK8_4,XK9_4,XK10_4,XK11_4,XK12_4,A1_1,A2_1,A1_2,A2_2,MU1,MU2,CU1,CU2,LAM1,LAM2,C11,C12,C21,C22,C31,C32,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K5_1,N2K6_1,N2K7_1,N2K8_1,N2K9_1,N2K10_1,N2K11_1,N2K12_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K5_2,N2K6_2,N2K7_2,N2K8_2,N2K9_2,N2K10_2,N2K11_2,N2K12_2,AL_V1,AL_V2,AL_V3,AL_V4,AL_V5,AL_V6,AL_V7,AL_V8,AL_V9,AL_V10,AL_V11,AL_V12,AL_D1,AL_D2,AL_D3,AL_D4,AL_D5,AL_D6,AL_D7,AL_D8,AL_D9,AL_D10,AL_D11,AL_D12,TT_V1,TT_V2,TT_V3,TT_V4,TT_V5,TT_V6,TT_V7,TT_V8,TT_V9,TT_V10,TT_V11,TT_V12,TT_D1,TT_D2,TT_D3,TT_D4,TT_D5,TT_D6,TT_D7,TT_D8,TT_D9,TT_D10,TT_D11,TT_D12,NR1,NR2,NR3,NR4,NR5,NR6,NR7,NR8,NR9,NR10,NR11,NR12,VE_V1,VE_V2,VE_V3,VE_V4,VE_V5,VE_V6,VE_V7,VE_V8,VE_V9,VE_V10,VE_V11,VE_V12,VE_D1,VE_D2,VE_D3,VE_D4,VE_D5,VE_D6,VE_D7,VE_D8,VE_D9,VE_D10,VE_D11,VE_D12,F0V1,F0V2,F0V3,F0V4,F0V5,F0V6,F0V7,F0V8,F0V9,F0V10,F0V11,F0V12,F0D1,F0D2,F0D3,F0D4,F0D5,F0D6,F0D7,F0D8,F0D9,F0D10,F0D11,F0D12)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    12-Nov-2019 13:46:24

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
t26 = n2k12_1.*t2;
t27 = n2k12_2.*t3;
t28 = 1.0./tt_v1;
t29 = 1.0./xk1_1;
t30 = tt_v1+ve_v1;
t31 = 1.0./t30;
t32 = 1.0./al_v1;
t33 = xk1_1.^t32;
t34 = t33.*tt_v1;
t35 = t4+t5+1.0;
t36 = t35.*ve_v1;
t37 = t34+t36;
t38 = 1.0./tt_v2;
t39 = 1.0./xk2_1;
t40 = tt_v2+ve_v2;
t41 = 1.0./t40;
t42 = 1.0./al_v2;
t43 = xk2_1.^t42;
t44 = t43.*tt_v2;
t45 = t6+t7+1.0;
t46 = t45.*ve_v2;
t47 = t44+t46;
t48 = 1.0./tt_v3;
t49 = 1.0./xk3_1;
t50 = tt_v3+ve_v3;
t51 = 1.0./t50;
t52 = 1.0./al_v3;
t53 = xk3_1.^t52;
t54 = t53.*tt_v3;
t55 = t8+t9+1.0;
t56 = t55.*ve_v3;
t57 = t54+t56;
t58 = 1.0./tt_v4;
t59 = 1.0./xk4_1;
t60 = tt_v4+ve_v4;
t61 = 1.0./t60;
t62 = 1.0./al_v4;
t63 = xk4_1.^t62;
t64 = t63.*tt_v4;
t65 = t10+t11+1.0;
t66 = t65.*ve_v4;
t67 = t64+t66;
t68 = 1.0./tt_v5;
t69 = 1.0./xk5_1;
t70 = tt_v5+ve_v5;
t71 = 1.0./t70;
t72 = 1.0./al_v5;
t73 = xk5_1.^t72;
t74 = t73.*tt_v5;
t75 = t12+t13+1.0;
t76 = t75.*ve_v5;
t77 = t74+t76;
t78 = 1.0./tt_v6;
t79 = 1.0./xk6_1;
t80 = tt_v6+ve_v6;
t81 = 1.0./t80;
t82 = 1.0./al_v6;
t83 = xk6_1.^t82;
t84 = t83.*tt_v6;
t85 = t14+t15+1.0;
t86 = t85.*ve_v6;
t87 = t84+t86;
t88 = 1.0./tt_v7;
t89 = 1.0./xk7_1;
t90 = tt_v7+ve_v7;
t91 = 1.0./t90;
t92 = 1.0./al_v7;
t93 = xk7_1.^t92;
t94 = t93.*tt_v7;
t95 = t16+t17+1.0;
t96 = t95.*ve_v7;
t97 = t94+t96;
t98 = 1.0./tt_v8;
t99 = 1.0./xk8_1;
t100 = tt_v8+ve_v8;
t101 = 1.0./t100;
t102 = 1.0./al_v8;
t103 = xk8_1.^t102;
t104 = t103.*tt_v8;
t105 = t18+t19+1.0;
t106 = t105.*ve_v8;
t107 = t104+t106;
t108 = 1.0./tt_v9;
t109 = 1.0./xk9_1;
t110 = tt_v9+ve_v9;
t111 = 1.0./t110;
t112 = 1.0./al_v9;
t113 = xk9_1.^t112;
t114 = t113.*tt_v9;
t115 = t20+t21+1.0;
t116 = t115.*ve_v9;
t117 = t114+t116;
t118 = 1.0./tt_v10;
t119 = 1.0./xk10_1;
t120 = tt_v10+ve_v10;
t121 = 1.0./t120;
t122 = 1.0./al_v10;
t123 = xk10_1.^t122;
t124 = t123.*tt_v10;
t125 = t22+t23+1.0;
t126 = t125.*ve_v10;
t127 = t124+t126;
t128 = 1.0./tt_v11;
t129 = 1.0./xk11_1;
t130 = tt_v11+ve_v11;
t131 = 1.0./t130;
t132 = 1.0./al_v11;
t133 = xk11_1.^t132;
t134 = t133.*tt_v11;
t135 = t24+t25+1.0;
t136 = t135.*ve_v11;
t137 = t134+t136;
t138 = 1.0./tt_v12;
t139 = 1.0./xk12_1;
t140 = tt_v12+ve_v12;
t141 = 1.0./t140;
t142 = 1.0./al_v12;
t143 = xk12_1.^t142;
t144 = t143.*tt_v12;
t145 = t26+t27+1.0;
t146 = t145.*ve_v12;
t147 = t144+t146;
t148 = 1.0./f0d1;
t149 = 1.0./f0d4;
t150 = 1.0./f0d7;
t151 = 1.0./f0d8;
t152 = 1.0./f0d10;
t153 = 1.0./f0d11;
t154 = 1.0./f0d9;
t155 = 1.0./f0d6;
t156 = 1.0./f0d5;
t157 = 1.0./f0d3;
t158 = 1.0./f0d2;
t159 = f0v1.*t31.*t37.*t148;
t160 = tt_d2+ve_d2;
t161 = 1.0./t160;
t162 = 1.0./al_d2;
t163 = xk2_3.^t162;
t164 = t163.*tt_d2;
t165 = tt_d4+ve_d4;
t166 = 1.0./t165;
t167 = 1.0./al_d4;
t168 = xk4_3.^t167;
t169 = t168.*tt_d4;
t170 = f0v4.*t61.*t67.*t149;
t171 = tt_d5+ve_d5;
t172 = 1.0./t171;
t173 = 1.0./al_d5;
t174 = xk5_3.^t173;
t175 = t174.*tt_d5;
t176 = tt_d6+ve_d6;
t177 = 1.0./t176;
t178 = f0v7.*t91.*t97.*t150;
t179 = tt_d8+ve_d8;
t180 = 1.0./t179;
t181 = 1.0./al_d8;
t182 = xk8_3.^t181;
t183 = t182.*tt_d8;
t184 = f0v8.*t101.*t107.*t151;
t185 = tt_d9+ve_d9;
t186 = 1.0./t185;
t187 = 1.0./al_d9;
t188 = xk9_3.^t187;
t189 = t188.*tt_d9;
t190 = f0v10.*t121.*t127.*t152;
t191 = tt_d11+ve_d11;
t192 = 1.0./t191;
t193 = 1.0./al_d11;
t194 = xk11_3.^t193;
t195 = t194.*tt_d11;
t196 = f0v11.*t131.*t137.*t153;
t197 = tt_d12+ve_d12;
t198 = 1.0./t197;
t199 = 1.0./al_d12;
t200 = xk12_3.^t199;
t201 = t200.*tt_d12;
t202 = t141.*t147.*ve_d12;
t203 = t201+t202;
t204 = f0d12.*t153.*t198.*t203;
t205 = t196+t204;
t206 = t205.*ve_d11;
t207 = t195+t206;
t208 = f0d11.*t152.*t192.*t207;
t209 = t190+t208;
t210 = t209.*ve_d10;
t211 = 1.0./al_d10;
t212 = xk10_3.^t211;
t213 = t212.*tt_d10;
t214 = t210+t213;
t215 = tt_d10+ve_d10;
t216 = 1.0./t215;
t217 = f0d10.*t154.*t214.*t216;
t218 = f0v9.*t111.*t117.*t154;
t219 = t217+t218;
t220 = t219.*ve_d9;
t221 = t189+t220;
t222 = f0d9.*t151.*t186.*t221;
t223 = t184+t222;
t224 = t223.*ve_d8;
t225 = t183+t224;
t226 = f0d8.*t150.*t180.*t225;
t227 = t178+t226;
t228 = t227.*ve_d7;
t229 = 1.0./al_d7;
t230 = xk7_3.^t229;
t231 = t230.*tt_d7;
t232 = t228+t231;
t233 = tt_d7+ve_d7;
t234 = 1.0./t233;
t235 = f0d7.*t155.*t232.*t234;
t236 = f0v6.*t81.*t87.*t155;
t237 = t235+t236;
t238 = t237.*ve_d6;
t239 = 1.0./al_d6;
t240 = xk6_3.^t239;
t241 = t240.*tt_d6;
t242 = t238+t241;
t243 = f0d6.*t156.*t177.*t242;
t244 = f0v5.*t71.*t77.*t156;
t245 = t243+t244;
t246 = t245.*ve_d5;
t247 = t175+t246;
t248 = f0d5.*t149.*t172.*t247;
t249 = t170+t248;
t250 = t249.*ve_d4;
t251 = t169+t250;
t252 = f0d4.*t157.*t166.*t251;
t253 = f0v3.*t51.*t57.*t157;
t254 = t252+t253;
t255 = t254.*ve_d3;
t256 = 1.0./al_d3;
t257 = xk3_3.^t256;
t258 = t257.*tt_d3;
t259 = t255+t258;
t260 = tt_d3+ve_d3;
t261 = 1.0./t260;
t262 = f0d3.*t158.*t259.*t261;
t263 = f0v2.*t41.*t47.*t158;
t264 = t262+t263;
t265 = t264.*ve_d2;
t266 = t164+t265;
t267 = f0d2.*t148.*t161.*t266;
t268 = 1.0./tt_d1;
t269 = tt_d1+ve_d1;
t270 = 1.0./t269;
t271 = 1.0./al_d1;
t272 = xk1_3.^t271;
t273 = t272.*tt_d1;
t274 = t159+t267;
t275 = t274.*ve_d1;
t276 = t273+t275;
t277 = 1.0./tt_d2;
t278 = 1.0./xk2_3;
t279 = 1.0./tt_d3;
t280 = 1.0./xk3_3;
t281 = 1.0./tt_d4;
t282 = 1.0./xk4_3;
t283 = 1.0./tt_d5;
t284 = 1.0./xk5_3;
t285 = 1.0./tt_d6;
t286 = 1.0./xk6_3;
t287 = 1.0./tt_d7;
t288 = 1.0./xk7_3;
t289 = 1.0./tt_d8;
t290 = 1.0./xk8_3;
t291 = 1.0./tt_d9;
t292 = 1.0./xk9_3;
t293 = 1.0./tt_d10;
t294 = 1.0./xk10_3;
t295 = 1.0./tt_d11;
t296 = 1.0./xk11_3;
t297 = 1.0./xk12_3;
t298 = 1.0./tt_d12;
f_s = [cu1+a1_1.*xn1_1+a1_2.*xn2_1-mu1.*xn1_2;cu2+a2_1.*xn1_1+a2_2.*xn2_1-mu2.*xn2_2;lam1.*(xn1_1-xn1_2);lam2.*(xn2_1-xn2_2);xn1_1-c11.*xn1_3;xn2_1-c12.*xn2_3;-(c31.*t2-c21.*xn1_3)./xn1_4;-(c32.*t3-c22.*xn2_3)./xn2_4;t28.*t29.*(t4+t5-t31.*t37+1.0);t38.*t39.*(t6+t7-t41.*t47+1.0);t48.*t49.*(t8+t9-t51.*t57+1.0);t58.*t59.*(t10+t11-t61.*t67+1.0);t68.*t69.*(t12+t13-t71.*t77+1.0);t78.*t79.*(t14+t15-t81.*t87+1.0);t88.*t89.*(t16+t17-t91.*t97+1.0);t98.*t99.*(t18+t19-t101.*t107+1.0);t108.*t109.*(t20+t21-t111.*t117+1.0);t118.*t119.*(t22+t23-t121.*t127+1.0);t128.*t129.*(t24+t25-t131.*t137+1.0);t138.*t139.*(t26+t27-t141.*t147+1.0);(t28.*((nr1+t4+t5)./nr1-t29.*t31.*t37.*xk1_2))./xk1_2;(t38.*((nr2+t6+t7)./nr2-t39.*t41.*t47.*xk2_2))./xk2_2;(t48.*((nr3+t8+t9)./nr3-t49.*t51.*t57.*xk3_2))./xk3_2;(t58.*((nr4+t10+t11)./nr4-t59.*t61.*t67.*xk4_2))./xk4_2;(t68.*((nr5+t12+t13)./nr5-t69.*t71.*t77.*xk5_2))./xk5_2;(t78.*((nr6+t14+t15)./nr6-t79.*t81.*t87.*xk6_2))./xk6_2;(t88.*((nr7+t16+t17)./nr7-t89.*t91.*t97.*xk7_2))./xk7_2;(t98.*((nr8+t18+t19)./nr8-t99.*t101.*t107.*xk8_2))./xk8_2;(t108.*((nr9+t20+t21)./nr9-t109.*t111.*t117.*xk9_2))./xk9_2;(t118.*((nr10+t22+t23)./nr10-t119.*t121.*t127.*xk10_2))./xk10_2;(t128.*((nr11+t24+t25)./nr11-t129.*t131.*t137.*xk11_2))./xk11_2;(t138.*((nr12+t26+t27)./nr12-t139.*t141.*t147.*xk12_2))./xk12_2;(t268.*(t159+t267-t270.*t276))./ve_d1;(t277.*(t262+t263-t161.*t266))./ve_d2;(t279.*(t252+t253-t259.*t261))./ve_d3;(t281.*(t170+t248-t166.*t251))./ve_d4;(t283.*(t243+t244-t172.*t247))./ve_d5;(t285.*(t235+t236-t177.*t242))./ve_d6;(t287.*(t178+t226-t232.*t234))./ve_d7;(t289.*(t184+t222-t180.*t225))./ve_d8;(t291.*(t217+t218-t186.*t221))./ve_d9;(t293.*(t190+t208-t214.*t216))./ve_d10;(t295.*(t196+t204-t192.*t207))./ve_d11;-t297.*t298.*(t198.*t203-(f0v12.*t141.*t147)./f0d12);(t268.*(-(t270.*t276.*xk1_4)./xk1_3+f0d2.*t148.*t161.*t266.*t278.*xk2_4+f0v1.*t29.*t31.*t37.*t148.*xk1_2))./xk1_4;(t277.*(-t161.*t266.*t278.*xk2_4+f0d3.*t158.*t259.*t261.*t280.*xk3_4+f0v2.*t39.*t41.*t47.*t158.*xk2_2))./xk2_4;(t279.*(-t259.*t261.*t280.*xk3_4+f0d4.*t157.*t166.*t251.*t282.*xk4_4+f0v3.*t49.*t51.*t57.*t157.*xk3_2))./xk3_4;(t281.*(-t166.*t251.*t282.*xk4_4+f0d5.*t149.*t172.*t247.*t284.*xk5_4+f0v4.*t59.*t61.*t67.*t149.*xk4_2))./xk4_4;(t283.*(-t172.*t247.*t284.*xk5_4+f0d6.*t156.*t177.*t242.*t286.*xk6_4+f0v5.*t69.*t71.*t77.*t156.*xk5_2))./xk5_4;(t285.*(-t177.*t242.*t286.*xk6_4+f0d7.*t155.*t232.*t234.*t288.*xk7_4+f0v6.*t79.*t81.*t87.*t155.*xk6_2))./xk6_4;(t287.*(-t232.*t234.*t288.*xk7_4+f0d8.*t150.*t180.*t225.*t290.*xk8_4+f0v7.*t89.*t91.*t97.*t150.*xk7_2))./xk7_4;(t289.*(-t180.*t225.*t290.*xk8_4+f0d9.*t151.*t186.*t221.*t292.*xk9_4+f0v8.*t99.*t101.*t107.*t151.*xk8_2))./xk8_4;(t291.*(-t186.*t221.*t292.*xk9_4+f0d10.*t154.*t214.*t216.*t294.*xk10_4+f0v9.*t109.*t111.*t117.*t154.*xk9_2))./xk9_4;(t293.*(-t214.*t216.*t294.*xk10_4+f0d11.*t152.*t192.*t207.*t296.*xk11_4+f0v10.*t119.*t121.*t127.*t152.*xk10_2))./xk10_4;(t295.*(-t192.*t207.*t296.*xk11_4+f0d12.*t153.*t198.*t203.*t297.*xk12_4+f0v11.*t129.*t131.*t137.*t153.*xk11_2))./xk11_4;(t298.*(t139.*t141.*t147.*xk12_2-t198.*t203.*t297.*xk12_4))./xk12_4];
