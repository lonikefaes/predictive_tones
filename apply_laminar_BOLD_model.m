% This script applies laminar BOLD model of Havlicek et al. (2021). Eight
% neuronal models are defined, each signifying target of modulatory input
% (superficial, middle or deep layers). The BOLD responses are simulated
% for multiple depths.

% Input: event related averages of BOLD data, depth estimates for each
% voxel. Stimulus is 4 tones (400 ms) separated by

% Works with single ROI for now

% Written by: Martin Havlicek, Isma Zulfiqar

clear;
close all;

Path = '\...\';     %Set own path to where scripts are located.
addpath(genpath(Path));


% ========= loading test dataset for a single ROI (Planum Polare, PP)
% normalized event related averages, at 0.8 isotropic; mask indices has
% corresponding voxel numbers in low-res (0.8), high-res (0.4) for the RoI
load('PP_LH_testData.mat'); 

% depth profile for each voxel, at 0.4 isotropic
% the depth estimates are made for a larger ROI, so no funny boundary
% effects happen for actual data points (smaller ROI)

load('PP_LH_layer_dist_testData.mat')
% ========= setting up model

condition_names = {'Comp_H', 'Comp_L', 'Oddball_H', 'Oddball_L', 'Unexp_H', 'Unexp_L'};
cmp1 = 1; % baseline
cmp2 = 3; % this condition

N = 3; % number of neuronal layers (usually 3), fixed
Ki = [7, 9]; % bold depths, can be varied. We have used [7, 9, 10, 11], takes long for each vascular depth so two as example.

TR = 1.6;

outfile = 'output_testData';


%%   loading ROI
Ana = T1.data;

close all
figure;
% ED and EV indicate depth estimates based on equidistant and equivolume
% approaches
% depth estimates for a single slice of RoI 
k = 60;  figure,overlayPlot(Ana(:,:,k),EVV2(:,:,k),1,0.3); axis image; colormap(gray); axis off;  title(num2str(k)); box off;
%k = 60;  figure,overlayPlot(Ana(:,:,k),EDV2(:,:,k),1,0.3); axis image; colormap(gray); axis off;  title(num2str(k)); box off;


EV_vox_r = EVV2(:);
%ED_vox_r = EDV2(:);

mat_ind = reshape([1:ceil(nx0/2)*ceil(ny0/2)*ceil(nz0/2)],ceil(nx0/2),ceil(ny0/2),ceil(nz0/2));
MAT_ind = [];

interp_factor = 2; % downsampling to get better estimate
for i = 1:ceil(nz0/2)
    MAT_ind = cat(3,MAT_ind,repmat(kron(mat_ind(:,:,i),ones(interp_factor)),[1,1,interp_factor]));
end

MAT_ind = MAT_ind(1:size(EVV2,1),1:size(EVV2,2),1:size(EVV2,3));

EV_vox_d = []; ED_vox_d = [];

% % downsampling...
% this takes time so might consider saving and then loading from file
disp('Downsampling depth estimates');
% for i = 1:length(unique(MAT_ind(:)))
%     list = MAT_ind==i;
%     if isempty(list)
%         disp(num2str(i));
%     end
%     EV_vox_d(i) = nanmean(EV_vox_r(MAT_ind==i));
%     %ED_vox_d(i) = nanmean(ED_vox_r(MAT_ind==i));
%     if rem(i,5000) == 0
%         fprintf('.');
%     end
% end
load('ds_vox_depth_PP_LH.mat')
disp(' Done!');

ind_vox = MAT_ind(mask_indices(:,3)); % getting mask indices in downsampled format
depth_map_HR = EVV2;
depth_map    = EV_vox_d;
sel = EV_vox_d(ind_vox);
sel_sel = sel>=0.0005 & sel<=0.9995;
ind_vox_sel = ind_vox(sel_sel);

depth_map_vox_sel = EV_vox_d(ind_vox_sel);


%% organising functional data (ER averages), per condition, all repetitions
% sound onset at 2nd timepoint (marked as 0)
perVoxResp = [];

figure(1000); clf
for cond = 1:6
    perVoxResp(cond,:,:) = eval(['squeeze(nanmean(nanmean(ER_avg_cond' num2str(cond) ',2),1));']);
    plot(squeeze(nanmean(perVoxResp(cond,:,:),2)), 'linewidth', 2)
    hold on;
end
xx = {-1:6};
x = gca; x.XTick = 1:8; x.XTickLabel = xx; x.XLim = [0 9]; xlabel('TR'); ylabel('signal change (%)');
legend(condition_names);
title('Avg response across ROI'); grid on;

gap_size = 20; %% we add these timepoints (some noise) between the two conditions, the TR is small
                % and we dont have enough timepoints between conditions, these are
                % disregarded during EM

%% Applying model
for d = 1:length(Ki)
    K = Ki(d);
    disp(['Running DCM for number of depths: ' num2str(K)]);
    
    layer_axr = linspace(0,1,2*(K)+1);
    layer_axr = layer_axr(2:2:end);
    
    % Sample data from voxel space to K-number of BOLD depths:
    
    figure(1001);
    for condition = 1:6
        eval(['mean_Data' num2str(condition) '= squeeze(perVoxResp(' num2str(condition) ',:,:));']);
        eval(['mean_Data' num2str(condition) '_vox_sel = mean_Data' num2str(condition) '(sel_sel,:);']);
        
        M.K  = K; % number of BOLD depths
        
        % Sample data from voxel space to K-number of BOLD depths:
        eval(['y' num2str(condition) ' = BOLD_voxels2layers_flipdata(mean_Data' num2str(condition) '_vox_sel,depth_map_vox_sel,M.K);']);
        
        eval(['y' num2str(condition) '(1,:) = 0;']);
        subplot(3,2,condition); eval(['plot(y' num2str(condition) ')']); title(condition_names(condition));
        xx = {-1:6};
        x = gca; x.XLim = [0 9]; x.XTick = 1:8; x.XTickLabel = xx; xlabel('TR')
        x.YLim = [-1 3]; ylabel('Signal change (%)')
    end
    
    % adding noise between responses, NOT fitted
    mid = []; for ii = 1:K, mid(:,ii) = 0.1*randn(gap_size,1);end
    
    
    %--------------------------------------------------------------------------
    % Specify neuronal and hemodynamic model for inversion:
    %--------------------------------------------------------------------------

    M.N  = N;
    dt = 0.05;
    
    % concatenate in time dimension
    Y.y = eval(['[y' num2str(cmp1) ';mid;y' num2str(cmp2) '];']);
    
    
    
    onset{1}    = [{[47.9]}, {[1.6 46.4]}];  % input onset specification in sec (modulatory goes always as first ones)
    duration{1} = [{[0.1]}, {[1.6 1.6]}];  % input durations (in sec)
    
    cnam   = {'c1m','c2d'};  % Input names
    
    
    mask = zeros(size(Y.y)); mask([2:8,30:36],:) = 1; % this is a temporal mask to ensure fitting only happens for datapoints
                                                       % and not for the gap
                                                          
    %Ki = K;
    % Loop over different number of BOLD cortical depths (later calculate Bayesian averange of neuronal parameters)
    %for d = 1:length(Ki)
    M.K  = Ki(d); % number of BOLD depths
    
    
    rnam   = strsplit(num2str(1:M.K));
    nr     = size(Y.y,2);  % depths
    ns     = size(Y.y,1);  % time-points
    
    cutoff = Inf;   % if some low pass filtering has to be done ... otherwise Inf for none
    
    DCM  = create_SPM_file_for_DCM2(Y,ns,TR,onset,duration,cnam,cutoff,rnam,round(TR/dt));
    
    DCM.Y.X0 = DCM.Y.X0(:,[2:end]);   % confounds matrix
    % DCM.Y
    % get model parameters M.P0
    M   = LBR_param_priors(M.N,M.K,zeros(N,N),zeros(N,N),zeros(N,2));
    
    M.P0.s = 0;
    M.P0.V0t = 3;
    M.P0.nr  = 3;
    M.P0.al_v  = 0.35;
    M.P0.w_v  = 0.5;
    %M.P0.n = 3;
    M.P0.nb = 0;
    M.P0.tau_v_same = 1;
    M.P0.tau_d_same = 1;
    
    %M.P0.tau_d_in = 30;
    %M.P0.tau_d_de = 30;
    pE  = spm_unvec(spm_vec(M.P0)*0,M.P0);
    pE.nb   = 0;
    
    pE.A   = zeros(N);
    pE.C   = [zeros(N,2),ones(N,1)/4];
    spC    = spm_unvec(spm_vec(pE)*0,pE);
    %pE.Bmu =
    
    M.pE = pE;
    M.pC = diag(spm_vec(spC));
    
    M.TR     = TR;
    M.delays = ones(1,nr)*(TR/2);
    M.m     = nr;
    M.n     = length(M.x(:));
    M.l     = nr;
    M.dt    = DCM.U.dt;
    M.ns    = ns;
    M.asl   = 0;
    M.IS    = 'spm_int_IT';
    
    %
    M.f       = @LBR_model_fx;
    M.g       = @LBR_model_gx;
    
    %
    M.Mask = reshape(find(mask),7*2,M.K); % temporal mask, 7 timepoints for both conditions

    % estimates voxel blurring for given number of BOLD depths
        % save this kernel in case re-running analysis, takes long!
    %M.kernel = BOLD_estimate_laminar_PSF_3D(M.N,M.K,MAT_ind,depth_map_HR,depth_map,ind_vox_sel);
    load("PP_LH_PSFkernel.mat");
    M.kernel = kernel;
    
    spC       = spm_unvec(spm_vec(M.P0)*0,M.P0);
    
    % Neuronal parameters
    % means:
    pE.A      = zeros(N);
    pE.C      = [zeros(N,1),ones(N,1)];
    pE.Bmu    = [zeros(3,1)];
    pE.Blam   = [0];
    %pE.B      = cat(3,zeros(N),zeros(N),eye(N)/5);
    pE.B      = cat(3,eye(N)*0);
    pE.mu = -0.8;
    pE.lam = 1.8;
    
    % variances:
    spC.C       = [zeros(N,1),ones(N,1)]*exp(0);
    spC.nsig    = exp(-4);
    spC.mu      = exp(0)*0; % ====> not fitting to reduce undershoot
    spC.sigma   = exp(-2);
    spC.lam     = exp(-2)*0; % ====> not fitting to reduce undershoot
    
    spC.B     = cat(3,diag([1 0 0])*exp(0.5));  % sigma targeting superficial depths
    %   spC.Bmu    = [ones(3,1)*exp(-1)];
    spC.Bmu    = [[0;0;0]*exp(-1)];
    
    spC.Blam   = [zeros(1,1)];
    
    
    
    % hemodynamic parameters:
    % means:
    pE.tau_d_de = 0;
    pE.tau_d_in = 0;
    pE.s_d      = 0;
    % variances:
    spC.tau_d_de  = exp(-3)*0;
    spC.tau_d_in  = exp(-3)*0;
    spC.al_d      = exp(-5);
    spC.s_d     = exp(-1);     % slope of draining vein
    %
    M.pE = pE;
    M.pC = diag(spm_vec(spC));
    
    DCM.M     = M;
    
    % ================== models to be tested =============== %
    DCM_M1    = DCM;
    %call model inversion:
    inv_DCM_M1   = LBR_model_inversion_dyn_mask(DCM_M1);
    
    DCMinv{1,d} = inv_DCM_M1;
    
    % Specify (only) what is different in the second model (neuronal part):
    DCM_M2    = DCM;
    M2        = M;
    spC2      = spC;
    spC2.B     = cat(3,diag([0 0 1])*exp(0.5));  % sigma targeting bottom depths
    
    M2.pC       = diag(spm_vec(spC2));
    DCM_M2.M    = M2;
    inv_DCM_M2   = LBR_model_inversion_dyn_mask(DCM_M2);
    DCMinv{2,d} = inv_DCM_M2;
    
    
    % Specify (only) what is different in the third model (neuronal part):
    DCM_M3    = DCM;
    M3        = M;
    spC3      = spC;
    spC3.B     = cat(3,diag([0 1 0])*exp(0.5));  % sigma targeting middle depths
    
    M3.pC       = diag(spm_vec(spC3));
    DCM_M3.M    = M3;
    inv_DCM_M3   = LBR_model_inversion_dyn_mask(DCM_M3);
    DCMinv{3,d} = inv_DCM_M3;
    
    % Specify (only) what is different in the fourth model (neuronal part):
    DCM_M4    = DCM;
    M4        = M;
    spC4      = spC;
    spC4.B     = cat(3,diag([1 1 1])*exp(0.5));  % sigma targeting all depths
    
    M4.pC       = diag(spm_vec(spC4));
    DCM_M4.M    = M4;
    inv_DCM_M4   = LBR_model_inversion_dyn_mask(DCM_M4);
    DCMinv{4,d} = inv_DCM_M4;
    
    % Specify (only) what is different in the fifth model (neuronal part):
    DCM_M5    = DCM;
    M5        = M;
    spC5      = spC;
    spC5.B     = cat(3,diag([1 0 1])*exp(0.5));  % sigma targeting top and bottom depths
    
    M5.pC       = diag(spm_vec(spC5));
    DCM_M5.M    = M5;
    inv_DCM_M5   = LBR_model_inversion_dyn_mask(DCM_M5);
    DCMinv{5,d} = inv_DCM_M5;
    
    
    % Specify (only) what is different in the sixth model (neuronal part):
    DCM_M6    = DCM;
    M6        = M;
    spC6      = spC;
    spC6.B     = cat(3,diag([1 1 0])*exp(0.5));  % sigma targeting top and mid depths
    
    M6.pC       = diag(spm_vec(spC6));
    DCM_M6.M    = M6;
    inv_DCM_M6   = LBR_model_inversion_dyn_mask(DCM_M6);
    DCMinv{6,d} = inv_DCM_M6;
    
    % Specify (only) what is different in the seventh model (neuronal part):
    DCM_M7    = DCM;
    M7        = M;
    spC7      = spC;
    spC7.B     = cat(3,diag([0 1 1])*exp(0.5));  % sigma targeting middle and bottom depths
    
    M7.pC       = diag(spm_vec(spC7));
    DCM_M7.M    = M7;
    inv_DCM_M7   = LBR_model_inversion_dyn_mask(DCM_M7);
    DCMinv{7,d} = inv_DCM_M7;
    
    % Specify (only) what is different in the eighth model (neuronal part):
    DCM_M8    = DCM;
    M8        = M;
    spC8      = spC;
    spC8.B     = cat(3,diag([0 0 0])*exp(0.5));  % no  modulation
    % add another level of B if additional mod condition
    M8.pC       = diag(spm_vec(spC8));
    DCM_M8.M    = M8;
    inv_DCM_M8   = LBR_model_inversion_dyn_mask(DCM_M8);
    DCMinv{8,d} = inv_DCM_M8;
    
    %============== changing Bmu, not using as no BOLD resolution in
    %time doesn't reflect Bmu
    %
    %     DCM_M9    = DCM;
    %     M9        = M;
    %     spC9      = spC;
    %     spC9.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC9.Bmu    = [[1;0;0]*exp(-1)]; % modulating top deths
    %
    %     M9.pC       = diag(spm_vec(spC9));
    %     DCM_M9.M    = M9;
    %     inv_DCM_M9   = LBR_model_inversion_dyn_mask(DCM_M9);
    %     DCMinv{9,d} = inv_DCM_M9;
    %
    %     % ====
    %
    %     DCM_M10    = DCM;
    %     M10        = M;
    %     spC10      = spC;
    %     spC10.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC10.Bmu    = [[0;1;0]*exp(-1)]; % modulating mid depths
    %
    %     M10.pC       = diag(spm_vec(spC10));
    %     DCM_M10.M    = M10;
    %     inv_DCM_M10   = LBR_model_inversion_dyn_mask(DCM_M10);
    %     DCMinv{10,d} = inv_DCM_M10;
    %
    %     % ====
    %
    %     DCM_M11    = DCM;
    %     M11        = M;
    %     spC11      = spC;
    %     spC11.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC11.Bmu    = [[0;0;1]*exp(-1)]; % modulating deep depths
    %
    %     M11.pC       = diag(spm_vec(spC11));
    %     DCM_M11.M    = M11;
    %     inv_DCM_M11   = LBR_model_inversion_dyn_mask(DCM_M11);
    %     DCMinv{11,d} = inv_DCM_M11;
    %
    %     % ====
    %
    %     DCM_M12    = DCM;
    %     M12        = M;
    %     spC12      = spC;
    %     spC12.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC12.Bmu    = [[1;1;1]*exp(-1)]; % modulating deep depths
    %
    %     M12.pC       = diag(spm_vec(spC12));
    %     DCM_M12.M    = M12;
    %     inv_DCM_M12   = LBR_model_inversion_dyn_mask(DCM_M12);
    %     DCMinv{12,d} = inv_DCM_M12;
    %
    %     % ====
    %
    %     DCM_M13    = DCM;
    %     M13        = M;
    %     spC13      = spC;
    %     spC13.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC13.Bmu    = [[1;0;1]*exp(-1)]; %
    %
    %     M13.pC       = diag(spm_vec(spC13));
    %     DCM_M13.M    = M13;
    %     inv_DCM_M13   = LBR_model_inversion_dyn_mask(DCM_M13);
    %     DCMinv{13,d} = inv_DCM_M13;
    %
    %     % ====
    %
    %     DCM_M14    = DCM;
    %     M14        = M;
    %     spC14      = spC;
    %     spC14.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC14.Bmu    = [[0;1;1]*exp(-1)]; % modulating deep depths
    %
    %     M14.pC       = diag(spm_vec(spC14));
    %     DCM_M14.M    = M14;
    %     inv_DCM_M14   = LBR_model_inversion_dyn_mask(DCM_M14);
    %     DCMinv{14,d} = inv_DCM_M14;
    %
    %     % ====
    %
    %     DCM_M15    = DCM;
    %     M15       = M;
    %     spC15      = spC;
    %     spC15.B     = cat(3,diag([0 0 0])*exp(0.5));  % setting B to zero
    %     spC15.Bmu    = [[1;1;0]*exp(-1)]; % modulating deep depths
    %
    %     M15.pC       = diag(spm_vec(spC15));
    %     DCM_M15.M    = M15;
    %     inv_DCM_M15   = LBR_model_inversion_dyn_mask(DCM_M15);
    %     DCMinv{15,d} = inv_DCM_M15;
    
end

% parameter averaging
yp = {};
Xp = {};
for m = 1:size(DCMinv,1)
    for bold_res = 1:size(DCMinv,2)
        DCM2avg{bold_res} =  DCMinv{m,bold_res};
    end
    
    BPA = spm_dcm_bpa(DCM2avg);
    
    [yp{m} Xp{m} dFx{m}] = spm_int_IT(BPA.Ep,BPA.M,BPA.U);
    
    F(m)    = BPA.F;  % accumulated over different number of BOLD depths
    Pp(m,:) = diag(BPA.Pp.B(:,:,1));
    Vp(m,:) = diag(BPA.Vp.B(:,:,1));
    Ep(m,:) = diag(BPA.Ep.B(:,:,1));
    
    Bu_est(m,:) = BPA.Ep.Bmu;
    
    disp(m)
end

% PERFORM MODEL COMPARISON...

dF = F - min(F);
p  = exp(dF)/sum(exp(dF));


figure, % model comparison results (in terms of posterior probability)
bar(p)

[~,bm] = max(p); % best model

% plotting data and fit from best model (averaged) (single time point,
% around peak)
figure(202); clf
colors = colormap(jet(Ki(1)));
for i = 1:Ki(1)
    subplot(1,2,1)
    plot(Y.y(2:8,i), ':o', 'color', colors(i,:), 'linewidth',1); hold on
    plot(yp{1,bm(1)}(2:8,i), 'x', 'color', colors(i,:),'linewidth',2);
    title(condition_names{cmp1}); 
    x= gca; x.YLim = [-0.1 2];x.XTick = 0:Ki(1)+1; xlabel('depth (csf - wm)'); ylabel('signal change (%)');
    grid on;
    title('condition 1')
    
    subplot(1,2,2);
    plot(Y.y(30:36,i), ':o', 'color', colors(i,:),'linewidth',1); hold on
    plot(yp{1,bm(1)}(30:36,i), 'x', 'color', colors(i,:),'linewidth',1);
    title(condition_names{cmp2});  
    x= gca; x.YLim = [-0.1 2]; x.XTick = 0:Ki(1)+1; xlabel('depth (csf - wm)'); ylabel('signal change (%)');
    grid on;
    title('condition 2')
end
x = gca; x.XTick = 1:Ki(1); xlabel('depth (csf - wm)'); ylabel('signal change (%)');
grid on
legend({'Data', 'Estimate'},'TextColor','k' )

[~,bm] = max(p)

disp(['Best model is Model: ' num2str(bm) ' with modulation targeting:'  newline num2str(Ep(bm,:)) ...
    newline 'sp mid and deep layers respectively.'])


