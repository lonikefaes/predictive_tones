function DCMinv = LBR_model_inversion_dyn(DCM)



[Ep,Cp,Eh,F] = spm_nlsi_GN_laminar(DCM.M,DCM.U,DCM.Y);

[yp X] = spm_int_IT(Ep,DCM.M,DCM.U);

DCMinv.Ep = Ep;
DCMinv.Cp = Cp;
DCMinv.U  = DCM.U;
DCMinv.Eh = Eh;
DCMinv.F  = F;
DCMinv.Y  = DCM.Y;
DCMinv.M  = DCM.M;
DCMinv.Yp  = yp;
DCMinv.Xp  = X;

return;
