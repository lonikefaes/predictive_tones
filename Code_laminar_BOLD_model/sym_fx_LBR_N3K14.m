function f_s = sym_fx_LBR_N3K14(xn1_1,xn2_1,xn3_1,xn1_2,xn2_2,xn3_2,xn1_3,xn2_3,xn3_3,xn1_4,xn2_4,xn3_4,xk1_1,xk2_1,xk3_1,xk4_1,xk5_1,xk6_1,xk7_1,xk8_1,xk9_1,xk10_1,xk11_1,xk12_1,xk13_1,xk14_1,xk1_2,xk2_2,xk3_2,xk4_2,xk5_2,xk6_2,xk7_2,xk8_2,xk9_2,xk10_2,xk11_2,xk12_2,xk13_2,xk14_2,xk1_3,xk2_3,xk3_3,xk4_3,xk5_3,xk6_3,xk7_3,xk8_3,xk9_3,xk10_3,xk11_3,xk12_3,xk13_3,xk14_3,xk1_4,xk2_4,xk3_4,xk4_4,xk5_4,xk6_4,xk7_4,xk8_4,xk9_4,xk10_4,xk11_4,xk12_4,xk13_4,xk14_4,a1_1,a2_1,a3_1,a1_2,a2_2,a3_2,a1_3,a2_3,a3_3,mu1,mu2,mu3,cu1,cu2,cu3,lam1,lam2,lam3,c11,c12,c13,c21,c22,c23,c31,c32,c33,n2k1_1,n2k2_1,n2k3_1,n2k4_1,n2k5_1,n2k6_1,n2k7_1,n2k8_1,n2k9_1,n2k10_1,n2k11_1,n2k12_1,n2k13_1,n2k14_1,n2k1_2,n2k2_2,n2k3_2,n2k4_2,n2k5_2,n2k6_2,n2k7_2,n2k8_2,n2k9_2,n2k10_2,n2k11_2,n2k12_2,n2k13_2,n2k14_2,n2k1_3,n2k2_3,n2k3_3,n2k4_3,n2k5_3,n2k6_3,n2k7_3,n2k8_3,n2k9_3,n2k10_3,n2k11_3,n2k12_3,n2k13_3,n2k14_3,al_v1,al_v2,al_v3,al_v4,al_v5,al_v6,al_v7,al_v8,al_v9,al_v10,al_v11,al_v12,al_v13,al_v14,al_d1,al_d2,al_d3,al_d4,al_d5,al_d6,al_d7,al_d8,al_d9,al_d10,al_d11,al_d12,al_d13,al_d14,tt_v1,tt_v2,tt_v3,tt_v4,tt_v5,tt_v6,tt_v7,tt_v8,tt_v9,tt_v10,tt_v11,tt_v12,tt_v13,tt_v14,tt_d1,tt_d2,tt_d3,tt_d4,tt_d5,tt_d6,tt_d7,tt_d8,tt_d9,tt_d10,tt_d11,tt_d12,tt_d13,tt_d14,nr1,nr2,nr3,nr4,nr5,nr6,nr7,nr8,nr9,nr10,nr11,nr12,nr13,nr14,ve_v1,ve_v2,ve_v3,ve_v4,ve_v5,ve_v6,ve_v7,ve_v8,ve_v9,ve_v10,ve_v11,ve_v12,ve_v13,ve_v14,ve_d1,ve_d2,ve_d3,ve_d4,ve_d5,ve_d6,ve_d7,ve_d8,ve_d9,ve_d10,ve_d11,ve_d12,ve_d13,ve_d14,f0v1,f0v2,f0v3,f0v4,f0v5,f0v6,f0v7,f0v8,f0v9,f0v10,f0v11,f0v12,f0v13,f0v14,f0d1,f0d2,f0d3,f0d4,f0d5,f0d6,f0d7,f0d8,f0d9,f0d10,f0d11,f0d12,f0d13,f0d14)
%SYM_FX_LBR_N3K14
%    F_S = SYM_FX_LBR_N3K14(XN1_1,XN2_1,XN3_1,XN1_2,XN2_2,XN3_2,XN1_3,XN2_3,XN3_3,XN1_4,XN2_4,XN3_4,XK1_1,XK2_1,XK3_1,XK4_1,XK5_1,XK6_1,XK7_1,XK8_1,XK9_1,XK10_1,XK11_1,XK12_1,XK13_1,XK14_1,XK1_2,XK2_2,XK3_2,XK4_2,XK5_2,XK6_2,XK7_2,XK8_2,XK9_2,XK10_2,XK11_2,XK12_2,XK13_2,XK14_2,XK1_3,XK2_3,XK3_3,XK4_3,XK5_3,XK6_3,XK7_3,XK8_3,XK9_3,XK10_3,XK11_3,XK12_3,XK13_3,XK14_3,XK1_4,XK2_4,XK3_4,XK4_4,XK5_4,XK6_4,XK7_4,XK8_4,XK9_4,XK10_4,XK11_4,XK12_4,XK13_4,XK14_4,A1_1,A2_1,A3_1,A1_2,A2_2,A3_2,A1_3,A2_3,A3_3,MU1,MU2,MU3,CU1,CU2,CU3,LAM1,LAM2,LAM3,C11,C12,C13,C21,C22,C23,C31,C32,C33,N2K1_1,N2K2_1,N2K3_1,N2K4_1,N2K5_1,N2K6_1,N2K7_1,N2K8_1,N2K9_1,N2K10_1,N2K11_1,N2K12_1,N2K13_1,N2K14_1,N2K1_2,N2K2_2,N2K3_2,N2K4_2,N2K5_2,N2K6_2,N2K7_2,N2K8_2,N2K9_2,N2K10_2,N2K11_2,N2K12_2,N2K13_2,N2K14_2,N2K1_3,N2K2_3,N2K3_3,N2K4_3,N2K5_3,N2K6_3,N2K7_3,N2K8_3,N2K9_3,N2K10_3,N2K11_3,N2K12_3,N2K13_3,N2K14_3,AL_V1,AL_V2,AL_V3,AL_V4,AL_V5,AL_V6,AL_V7,AL_V8,AL_V9,AL_V10,AL_V11,AL_V12,AL_V13,AL_V14,AL_D1,AL_D2,AL_D3,AL_D4,AL_D5,AL_D6,AL_D7,AL_D8,AL_D9,AL_D10,AL_D11,AL_D12,AL_D13,AL_D14,TT_V1,TT_V2,TT_V3,TT_V4,TT_V5,TT_V6,TT_V7,TT_V8,TT_V9,TT_V10,TT_V11,TT_V12,TT_V13,TT_V14,TT_D1,TT_D2,TT_D3,TT_D4,TT_D5,TT_D6,TT_D7,TT_D8,TT_D9,TT_D10,TT_D11,TT_D12,TT_D13,TT_D14,NR1,NR2,NR3,NR4,NR5,NR6,NR7,NR8,NR9,NR10,NR11,NR12,NR13,NR14,VE_V1,VE_V2,VE_V3,VE_V4,VE_V5,VE_V6,VE_V7,VE_V8,VE_V9,VE_V10,VE_V11,VE_V12,VE_V13,VE_V14,VE_D1,VE_D2,VE_D3,VE_D4,VE_D5,VE_D6,VE_D7,VE_D8,VE_D9,VE_D10,VE_D11,VE_D12,VE_D13,VE_D14,F0V1,F0V2,F0V3,F0V4,F0V5,F0V6,F0V7,F0V8,F0V9,F0V10,F0V11,F0V12,F0V13,F0V14,F0D1,F0D2,F0D3,F0D4,F0D5,F0D6,F0D7,F0D8,F0D9,F0D10,F0D11,F0D12,F0D13,F0D14)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    22-Jan-2020 02:05:48

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
t26 = n2k8_1.*t2;
t27 = n2k8_2.*t3;
t28 = n2k8_3.*t4;
t29 = n2k9_1.*t2;
t30 = n2k9_2.*t3;
t31 = n2k9_3.*t4;
t32 = n2k10_1.*t2;
t33 = n2k10_2.*t3;
t34 = n2k10_3.*t4;
t35 = n2k11_1.*t2;
t36 = n2k11_2.*t3;
t37 = n2k11_3.*t4;
t38 = n2k12_1.*t2;
t39 = n2k12_2.*t3;
t40 = n2k12_3.*t4;
t41 = n2k13_1.*t2;
t42 = n2k13_2.*t3;
t43 = n2k13_3.*t4;
t44 = n2k14_1.*t2;
t45 = n2k14_2.*t3;
t46 = n2k14_3.*t4;
t47 = 1.0./tt_v1;
t48 = 1.0./xk1_1;
t49 = tt_v1+ve_v1;
t50 = 1.0./t49;
t51 = t5+t6+t7+1.0;
t52 = t51.*ve_v1;
t53 = 1.0./al_v1;
t54 = xk1_1.^t53;
t55 = t54.*tt_v1;
t56 = t52+t55;
t57 = 1.0./tt_v2;
t58 = 1.0./xk2_1;
t59 = tt_v2+ve_v2;
t60 = 1.0./t59;
t61 = t8+t9+t10+1.0;
t62 = t61.*ve_v2;
t63 = 1.0./al_v2;
t64 = xk2_1.^t63;
t65 = t64.*tt_v2;
t66 = t62+t65;
t67 = 1.0./tt_v3;
t68 = 1.0./xk3_1;
t69 = tt_v3+ve_v3;
t70 = 1.0./t69;
t71 = t11+t12+t13+1.0;
t72 = t71.*ve_v3;
t73 = 1.0./al_v3;
t74 = xk3_1.^t73;
t75 = t74.*tt_v3;
t76 = t72+t75;
t77 = 1.0./tt_v4;
t78 = 1.0./xk4_1;
t79 = tt_v4+ve_v4;
t80 = 1.0./t79;
t81 = t14+t15+t16+1.0;
t82 = t81.*ve_v4;
t83 = 1.0./al_v4;
t84 = xk4_1.^t83;
t85 = t84.*tt_v4;
t86 = t82+t85;
t87 = 1.0./tt_v5;
t88 = 1.0./xk5_1;
t89 = tt_v5+ve_v5;
t90 = 1.0./t89;
t91 = t17+t18+t19+1.0;
t92 = t91.*ve_v5;
t93 = 1.0./al_v5;
t94 = xk5_1.^t93;
t95 = t94.*tt_v5;
t96 = t92+t95;
t97 = 1.0./tt_v6;
t98 = 1.0./xk6_1;
t99 = tt_v6+ve_v6;
t100 = 1.0./t99;
t101 = t20+t21+t22+1.0;
t102 = t101.*ve_v6;
t103 = 1.0./al_v6;
t104 = xk6_1.^t103;
t105 = t104.*tt_v6;
t106 = t102+t105;
t107 = 1.0./tt_v7;
t108 = 1.0./xk7_1;
t109 = tt_v7+ve_v7;
t110 = 1.0./t109;
t111 = t23+t24+t25+1.0;
t112 = t111.*ve_v7;
t113 = 1.0./al_v7;
t114 = xk7_1.^t113;
t115 = t114.*tt_v7;
t116 = t112+t115;
t117 = 1.0./tt_v8;
t118 = 1.0./xk8_1;
t119 = tt_v8+ve_v8;
t120 = 1.0./t119;
t121 = t26+t27+t28+1.0;
t122 = t121.*ve_v8;
t123 = 1.0./al_v8;
t124 = xk8_1.^t123;
t125 = t124.*tt_v8;
t126 = t122+t125;
t127 = 1.0./tt_v9;
t128 = 1.0./xk9_1;
t129 = tt_v9+ve_v9;
t130 = 1.0./t129;
t131 = t29+t30+t31+1.0;
t132 = t131.*ve_v9;
t133 = 1.0./al_v9;
t134 = xk9_1.^t133;
t135 = t134.*tt_v9;
t136 = t132+t135;
t137 = 1.0./tt_v10;
t138 = 1.0./xk10_1;
t139 = tt_v10+ve_v10;
t140 = 1.0./t139;
t141 = t32+t33+t34+1.0;
t142 = t141.*ve_v10;
t143 = 1.0./al_v10;
t144 = xk10_1.^t143;
t145 = t144.*tt_v10;
t146 = t142+t145;
t147 = 1.0./tt_v11;
t148 = 1.0./xk11_1;
t149 = tt_v11+ve_v11;
t150 = 1.0./t149;
t151 = t35+t36+t37+1.0;
t152 = t151.*ve_v11;
t153 = 1.0./al_v11;
t154 = xk11_1.^t153;
t155 = t154.*tt_v11;
t156 = t152+t155;
t157 = 1.0./tt_v12;
t158 = 1.0./xk12_1;
t159 = tt_v12+ve_v12;
t160 = 1.0./t159;
t161 = t38+t39+t40+1.0;
t162 = t161.*ve_v12;
t163 = 1.0./al_v12;
t164 = xk12_1.^t163;
t165 = t164.*tt_v12;
t166 = t162+t165;
t167 = 1.0./tt_v13;
t168 = 1.0./xk13_1;
t169 = tt_v13+ve_v13;
t170 = 1.0./t169;
t171 = t41+t42+t43+1.0;
t172 = t171.*ve_v13;
t173 = 1.0./al_v13;
t174 = xk13_1.^t173;
t175 = t174.*tt_v13;
t176 = t172+t175;
t177 = 1.0./tt_v14;
t178 = 1.0./xk14_1;
t179 = tt_v14+ve_v14;
t180 = 1.0./t179;
t181 = t44+t45+t46+1.0;
t182 = t181.*ve_v14;
t183 = 1.0./al_v14;
t184 = xk14_1.^t183;
t185 = t184.*tt_v14;
t186 = t182+t185;
t187 = 1.0./f0d1;
t188 = 1.0./f0d2;
t189 = 1.0./f0d4;
t190 = 1.0./f0d5;
t191 = 1.0./f0d6;
t192 = 1.0./f0d7;
t193 = 1.0./f0d8;
t194 = 1.0./f0d9;
t195 = 1.0./f0d10;
t196 = 1.0./f0d11;
t197 = 1.0./f0d12;
t198 = 1.0./f0d13;
t199 = 1.0./f0d3;
t200 = f0v1.*t50.*t56.*t187;
t201 = tt_d2+ve_d2;
t202 = 1.0./t201;
t203 = 1.0./al_d2;
t204 = xk2_3.^t203;
t205 = t204.*tt_d2;
t206 = f0v2.*t60.*t66.*t188;
t207 = tt_d4+ve_d4;
t208 = 1.0./t207;
t209 = f0v4.*t80.*t86.*t189;
t210 = tt_d5+ve_d5;
t211 = 1.0./t210;
t212 = 1.0./al_d5;
t213 = xk5_3.^t212;
t214 = t213.*tt_d5;
t215 = f0v5.*t90.*t96.*t190;
t216 = f0v6.*t100.*t106.*t191;
t217 = tt_d7+ve_d7;
t218 = 1.0./t217;
t219 = 1.0./al_d7;
t220 = xk7_3.^t219;
t221 = t220.*tt_d7;
t222 = f0v7.*t110.*t116.*t192;
t223 = tt_d8+ve_d8;
t224 = 1.0./t223;
t225 = 1.0./al_d8;
t226 = xk8_3.^t225;
t227 = t226.*tt_d8;
t228 = f0v8.*t120.*t126.*t193;
t229 = f0v9.*t130.*t136.*t194;
t230 = tt_d10+ve_d10;
t231 = 1.0./t230;
t232 = 1.0./al_d10;
t233 = xk10_3.^t232;
t234 = t233.*tt_d10;
t235 = f0v10.*t140.*t146.*t195;
t236 = tt_d11+ve_d11;
t237 = 1.0./t236;
t238 = 1.0./al_d11;
t239 = xk11_3.^t238;
t240 = t239.*tt_d11;
t241 = f0v11.*t150.*t156.*t196;
t242 = tt_d12+ve_d12;
t243 = 1.0./t242;
t244 = f0v12.*t160.*t166.*t197;
t245 = tt_d13+ve_d13;
t246 = 1.0./t245;
t247 = 1.0./al_d13;
t248 = xk13_3.^t247;
t249 = t248.*tt_d13;
t250 = f0v13.*t170.*t176.*t198;
t251 = 1.0./al_d14;
t252 = xk14_3.^t251;
t253 = t252.*tt_d14;
t254 = t180.*t186.*ve_d14;
t255 = t253+t254;
t256 = tt_d14+ve_d14;
t257 = 1.0./t256;
t258 = f0d14.*t198.*t255.*t257;
t259 = t250+t258;
t260 = t259.*ve_d13;
t261 = t249+t260;
t262 = f0d13.*t197.*t246.*t261;
t263 = t244+t262;
t264 = t263.*ve_d12;
t265 = 1.0./al_d12;
t266 = xk12_3.^t265;
t267 = t266.*tt_d12;
t268 = t264+t267;
t269 = f0d12.*t196.*t243.*t268;
t270 = t241+t269;
t271 = t270.*ve_d11;
t272 = t240+t271;
t273 = f0d11.*t195.*t237.*t272;
t274 = t235+t273;
t275 = t274.*ve_d10;
t276 = t234+t275;
t277 = f0d10.*t194.*t231.*t276;
t278 = t229+t277;
t279 = t278.*ve_d9;
t280 = 1.0./al_d9;
t281 = xk9_3.^t280;
t282 = t281.*tt_d9;
t283 = t279+t282;
t284 = tt_d9+ve_d9;
t285 = 1.0./t284;
t286 = f0d9.*t193.*t283.*t285;
t287 = t228+t286;
t288 = t287.*ve_d8;
t289 = t227+t288;
t290 = f0d8.*t192.*t224.*t289;
t291 = t222+t290;
t292 = t291.*ve_d7;
t293 = t221+t292;
t294 = f0d7.*t191.*t218.*t293;
t295 = t216+t294;
t296 = t295.*ve_d6;
t297 = 1.0./al_d6;
t298 = xk6_3.^t297;
t299 = t298.*tt_d6;
t300 = t296+t299;
t301 = tt_d6+ve_d6;
t302 = 1.0./t301;
t303 = f0d6.*t190.*t300.*t302;
t304 = t215+t303;
t305 = t304.*ve_d5;
t306 = t214+t305;
t307 = f0d5.*t189.*t211.*t306;
t308 = t209+t307;
t309 = t308.*ve_d4;
t310 = 1.0./al_d4;
t311 = xk4_3.^t310;
t312 = t311.*tt_d4;
t313 = t309+t312;
t314 = f0d4.*t199.*t208.*t313;
t315 = f0v3.*t70.*t76.*t199;
t316 = t314+t315;
t317 = t316.*ve_d3;
t318 = 1.0./al_d3;
t319 = xk3_3.^t318;
t320 = t319.*tt_d3;
t321 = t317+t320;
t322 = tt_d3+ve_d3;
t323 = 1.0./t322;
t324 = f0d3.*t188.*t321.*t323;
t325 = t206+t324;
t326 = t325.*ve_d2;
t327 = t205+t326;
t328 = f0d2.*t187.*t202.*t327;
t329 = 1.0./tt_d1;
t330 = tt_d1+ve_d1;
t331 = 1.0./t330;
t332 = t200+t328;
t333 = t332.*ve_d1;
t334 = 1.0./al_d1;
t335 = xk1_3.^t334;
t336 = t335.*tt_d1;
t337 = t333+t336;
t338 = 1.0./tt_d2;
t339 = 1.0./xk2_3;
t340 = 1.0./tt_d3;
t341 = 1.0./xk3_3;
t342 = 1.0./tt_d4;
t343 = 1.0./xk4_3;
t344 = 1.0./tt_d5;
t345 = 1.0./xk5_3;
t346 = 1.0./tt_d6;
t347 = 1.0./xk6_3;
t348 = 1.0./tt_d7;
t349 = 1.0./xk7_3;
t350 = 1.0./tt_d8;
t351 = 1.0./xk8_3;
t352 = 1.0./tt_d9;
t353 = 1.0./xk9_3;
t354 = 1.0./tt_d10;
t355 = 1.0./xk10_3;
t356 = 1.0./tt_d11;
t357 = 1.0./xk11_3;
t358 = 1.0./tt_d12;
t359 = 1.0./xk12_3;
t360 = 1.0./tt_d13;
t361 = 1.0./xk13_3;
t362 = 1.0./xk14_3;
t363 = 1.0./tt_d14;
f_s = [cu1+a1_1.*xn1_1+a1_2.*xn2_1+a1_3.*xn3_1-mu1.*xn1_2;cu2+a2_1.*xn1_1+a2_2.*xn2_1+a2_3.*xn3_1-mu2.*xn2_2;cu3+a3_1.*xn1_1+a3_2.*xn2_1+a3_3.*xn3_1-mu3.*xn3_2;lam1.*(xn1_1-xn1_2);lam2.*(xn2_1-xn2_2);lam3.*(xn3_1-xn3_2);xn1_1-c11.*xn1_3;xn2_1-c12.*xn2_3;xn3_1-c13.*xn3_3;-(c31.*t2-c21.*xn1_3)./xn1_4;-(c32.*t3-c22.*xn2_3)./xn2_4;-(c33.*t4-c23.*xn3_3)./xn3_4;t47.*t48.*(t5+t6+t7-t50.*t56+1.0);t57.*t58.*(t8+t9+t10-t60.*t66+1.0);t67.*t68.*(t11+t12+t13-t70.*t76+1.0);t77.*t78.*(t14+t15+t16-t80.*t86+1.0);t87.*t88.*(t17+t18+t19-t90.*t96+1.0);t97.*t98.*(t20+t21+t22-t100.*t106+1.0);t107.*t108.*(t23+t24+t25-t110.*t116+1.0);t117.*t118.*(t26+t27+t28-t120.*t126+1.0);t127.*t128.*(t29+t30+t31-t130.*t136+1.0);t137.*t138.*(t32+t33+t34-t140.*t146+1.0);t147.*t148.*(t35+t36+t37-t150.*t156+1.0);t157.*t158.*(t38+t39+t40-t160.*t166+1.0);t167.*t168.*(t41+t42+t43-t170.*t176+1.0);t177.*t178.*(t44+t45+t46-t180.*t186+1.0);(t47.*((nr1+t5+t6+t7)./nr1-t48.*t50.*t56.*xk1_2))./xk1_2;(t57.*((nr2+t8+t9+t10)./nr2-t58.*t60.*t66.*xk2_2))./xk2_2;(t67.*((nr3+t11+t12+t13)./nr3-t68.*t70.*t76.*xk3_2))./xk3_2;(t77.*((nr4+t14+t15+t16)./nr4-t78.*t80.*t86.*xk4_2))./xk4_2;(t87.*((nr5+t17+t18+t19)./nr5-t88.*t90.*t96.*xk5_2))./xk5_2;(t97.*((nr6+t20+t21+t22)./nr6-t98.*t100.*t106.*xk6_2))./xk6_2;(t107.*((nr7+t23+t24+t25)./nr7-t108.*t110.*t116.*xk7_2))./xk7_2;(t117.*((nr8+t26+t27+t28)./nr8-t118.*t120.*t126.*xk8_2))./xk8_2;(t127.*((nr9+t29+t30+t31)./nr9-t128.*t130.*t136.*xk9_2))./xk9_2;(t137.*((nr10+t32+t33+t34)./nr10-t138.*t140.*t146.*xk10_2))./xk10_2;(t147.*((nr11+t35+t36+t37)./nr11-t148.*t150.*t156.*xk11_2))./xk11_2;(t157.*((nr12+t38+t39+t40)./nr12-t158.*t160.*t166.*xk12_2))./xk12_2;(t167.*((nr13+t41+t42+t43)./nr13-t168.*t170.*t176.*xk13_2))./xk13_2;(t177.*((nr14+t44+t45+t46)./nr14-t178.*t180.*t186.*xk14_2))./xk14_2;(t329.*(t200+t328-t331.*t337))./ve_d1;(t338.*(t206+t324-t202.*t327))./ve_d2;(t340.*(t314+t315-t321.*t323))./ve_d3;(t342.*(t209+t307-t208.*t313))./ve_d4;(t344.*(t215+t303-t211.*t306))./ve_d5;(t346.*(t216+t294-t300.*t302))./ve_d6;(t348.*(t222+t290-t218.*t293))./ve_d7;(t350.*(t228+t286-t224.*t289))./ve_d8;(t352.*(t229+t277-t283.*t285))./ve_d9;(t354.*(t235+t273-t231.*t276))./ve_d10;(t356.*(t241+t269-t237.*t272))./ve_d11;(t358.*(t244+t262-t243.*t268))./ve_d12;(t360.*(t250+t258-t246.*t261))./ve_d13;-t362.*t363.*(t255.*t257-(f0v14.*t180.*t186)./f0d14);(t329.*(-(t331.*t337.*xk1_4)./xk1_3+f0d2.*t187.*t202.*t327.*t339.*xk2_4+f0v1.*t48.*t50.*t56.*t187.*xk1_2))./xk1_4;(t338.*(-t202.*t327.*t339.*xk2_4+f0d3.*t188.*t321.*t323.*t341.*xk3_4+f0v2.*t58.*t60.*t66.*t188.*xk2_2))./xk2_4;(t340.*(-t321.*t323.*t341.*xk3_4+f0d4.*t199.*t208.*t313.*t343.*xk4_4+f0v3.*t68.*t70.*t76.*t199.*xk3_2))./xk3_4;(t342.*(-t208.*t313.*t343.*xk4_4+f0d5.*t189.*t211.*t306.*t345.*xk5_4+f0v4.*t78.*t80.*t86.*t189.*xk4_2))./xk4_4;(t344.*(-t211.*t306.*t345.*xk5_4+f0d6.*t190.*t300.*t302.*t347.*xk6_4+f0v5.*t88.*t90.*t96.*t190.*xk5_2))./xk5_4;(t346.*(-t300.*t302.*t347.*xk6_4+f0d7.*t191.*t218.*t293.*t349.*xk7_4+f0v6.*t98.*t100.*t106.*t191.*xk6_2))./xk6_4;(t348.*(-t218.*t293.*t349.*xk7_4+f0d8.*t192.*t224.*t289.*t351.*xk8_4+f0v7.*t108.*t110.*t116.*t192.*xk7_2))./xk7_4;(t350.*(-t224.*t289.*t351.*xk8_4+f0d9.*t193.*t283.*t285.*t353.*xk9_4+f0v8.*t118.*t120.*t126.*t193.*xk8_2))./xk8_4;(t352.*(-t283.*t285.*t353.*xk9_4+f0d10.*t194.*t231.*t276.*t355.*xk10_4+f0v9.*t128.*t130.*t136.*t194.*xk9_2))./xk9_4;(t354.*(-t231.*t276.*t355.*xk10_4+f0d11.*t195.*t237.*t272.*t357.*xk11_4+f0v10.*t138.*t140.*t146.*t195.*xk10_2))./xk10_4;(t356.*(-t237.*t272.*t357.*xk11_4+f0d12.*t196.*t243.*t268.*t359.*xk12_4+f0v11.*t148.*t150.*t156.*t196.*xk11_2))./xk11_4;(t358.*(-t243.*t268.*t359.*xk12_4+f0d13.*t197.*t246.*t261.*t361.*xk13_4+f0v12.*t158.*t160.*t166.*t197.*xk12_2))./xk12_4;(t360.*(-t246.*t261.*t361.*xk13_4+f0d14.*t198.*t255.*t257.*t362.*xk14_4+f0v13.*t168.*t170.*t176.*t198.*xk13_2))./xk13_4;(t363.*(t178.*t180.*t186.*xk14_2-t255.*t257.*t362.*xk14_4))./xk14_4];
