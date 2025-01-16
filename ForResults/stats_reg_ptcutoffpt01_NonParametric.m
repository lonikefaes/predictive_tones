%%% Non-parametric testing, averaged over hemisphere and cond

%% ===================================================================== %%
%%  reg, t = 2 , unflipped, separate rois, data cleaned hard cutoff pt01     %%
%% ===================================================================== %%

clear; close all;
%close all;

Path = '\PathtoData\';  %Insert path

subjects = {'S1','S2','S3','S5','S6','S7','S8','S9','S10','S11'};

rois = {'HG', 'PP', 'PT', 'aSTG', 'pSTG'};
hemispheres = {'LH','RH'};

conditions = {'H', 'L'};

nordic = 'reg'; % reg or nordic or NoN

vox_tag = '_ptcutoffpt01_notflipped';
depths = '_depths791011';


for s = 1:length(subjects)
    
    for r = 1:length(rois)
        
        for h = 1:length(hemispheres)
            
            inpath = [Path 'Data_from_Lonike\Output\' subjects{s} '\DynamicDCM\' rois{r} '_' hemispheres{h} '\results_Comp_'];
            
            for c = 1:length(conditions)
                %oddball
                load([inpath conditions{c} '_Oddball_' conditions{c} depths '_mask_' nordic vox_tag '.mat'], 'p', 'Ep', 'mask_indices');
                
                if length(mask_indices)>300
                    [x y] = max(p);
                    eval(['BM_oddball' conditions{c} '_' rois{r} '_' hemispheres{h} '(s,1:3) = round(Ep(y,:),2);']);
                else
                    disp([subjects{s} '_' rois{r} '_' hemispheres{h}]);
                    eval(['BM_oddball' conditions{c} '_' rois{r} '_' hemispheres{h} '(s,1:3) = nan(1,3);']);
                end
            end
        end
    end
end

%% Create datamatrix
datamat = []; 
for r = 1:length(rois)
    
    for h = 1:length(hemispheres)
        
        for c = 1:length(conditions)
            
            temp = eval(['BM_oddball' conditions{c} '_' rois{r} '_' hemispheres{h}]);
            temp = fliplr(temp); % deep mid sup
            
            datamat(:,:,r,h,c) = temp;
            
        end
    end
end

Permutations = zeros(2^10-2,10); 
count1       = 1;
for ind = 1:9
    % cycling from 1 to n-1 rather than 0 to n
    C = nchoosek(1:10,ind);
    P = ones(size(C,1),10);
    for indP = 1:size(C,1)
        P(indP,C(indP,:))=-1;
    end
    Permutations(count1 + (0:size(P,1)-1),:) = P;
    count1  = count1 + size(P,1);
end
Permutations(end+1,:) = repmat(-1,1,10); 
Permutations(end+1,:) = repmat(1,1,10);


%%  ==================== HG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,1,:,:)); % HG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); % We have added the trivial permutations, two-sided test.

pvaluesCorr(1,:) = p;
tval(1,:) = t_data;

resp_HG_oddball = datamat11;


%%  ==================== PP
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,2,:,:)); % PP

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(2,:) = p;
tval(2,:) = t_data;


resp_PP_oddball = datamat11;


%%  ==================== PT
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,3,:,:)); % PT

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat;  p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end
p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(3,:) = p;
tval(3,:) = t_data;

resp_PT_oddball = datamat11;

%%  ==================== aSTG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,4,:,:)); % aSTG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end
p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(4,:) = p;
tval(4,:) = t_data;

resp_aSTG_oddball = datamat11;

%%  ==================== pSTG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,5,:,:)); % pSTG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(5,:) = p;
tval(5,:) = t_data;

resp_pSTG_oddball = datamat11;

%% Calculate FDR

pval_vec = [pvaluesCorr(1,:), pvaluesCorr(2,:), pvaluesCorr(3,:), pvaluesCorr(4,:), pvaluesCorr(5,:)];

[FDR] = mafdr(pval_vec, 'BHFDR', 'true')       %Take all 15 numbers and correct.

pvaluesCorr = [FDR(1:3); FDR(4:6); FDR(7:9); FDR(10:12); FDR(13:15)];
pvaluesCorr = table(pvaluesCorr);
tvalTab = table(tval);

%writetable(pvaluesCorr, ['pvalCorr', depths, '_permutations_FDRCorr_Mispred.txt'], 'Delimiter', '\t')
%writetable(tvalTab, ['tval', depths, '_permutations_Mispred.txt'], 'Delimiter', '\t')

save(['output_t2_model_Mispred', vox_tag, depths, '.mat'], 'resp_HG_oddball', 'resp_PP_oddball', 'resp_PT_oddball', ...
   'resp_aSTG_oddball', 'resp_pSTG_oddball');

%% UNPREDICTABLE %%

%% ===================================================================== %%
%%  reg, t = 2 , unflipped, separate rois, data cleaned hardcutoff pt01     %%
%% ===================================================================== %%

clear; close all;
%close all;

Path = 'F:\Tones\Laminar_BOLD_model\';

%subjects = {'S1','S5','S6','S7','S8','S9','S10','S11'};
subjects = {'S1','S2', 'S3','S5','S6','S7','S8','S9','S10','S11'};

rois = {'HG', 'PP', 'PT', 'aSTG', 'pSTG'};
%rois = {'HG', 'PT', 'pSTG'};
hemispheres = {'LH','RH'};

conditions = {'H', 'L'};

nordic = 'reg'; % reg or nordic or NoN

vox_tag = '_ptcutoffpt01_notflipped';
depths = '_depths791011';


for s = 1:length(subjects)
    
    for r = 1:length(rois)
        
        for h = 1:length(hemispheres)
            
            inpath = [Path 'Data_from_Lonike\Output\' subjects{s} '\DynamicDCM\' rois{r} '_' hemispheres{h} '\results_Comp_'];
            
            for c = 1:length(conditions)
                %oddball
                load([inpath conditions{c} '_Unexp_' conditions{c} depths '_mask_' nordic vox_tag '.mat'], 'p', 'Ep', 'mask_indices');
                
                if length(mask_indices)>300
                    [x y] = max(p);
                    eval(['BM_unexp' conditions{c} '_' rois{r} '_' hemispheres{h} '(s,1:3) = round(Ep(y,:),2);']);
                else
                    disp([subjects{s} '_' rois{r} '_' hemispheres{h}]);
                    eval(['BM_unexp' conditions{c} '_' rois{r} '_' hemispheres{h} '(s,1:3) = nan(1,3);']);
                end
            end
        end
    end
end

%% Create datamatrix
datamat = []; 
for r = 1:length(rois)
    
    for h = 1:length(hemispheres)
        
        for c = 1:length(conditions)
            
            temp = eval(['BM_unexp' conditions{c} '_' rois{r} '_' hemispheres{h}]);
            temp = fliplr(temp); % deep mid sup
            
            datamat(:,:,r,h,c) = temp;
            
        end
    end
end

Permutations = zeros(2^10-2,10); 
count1       = 1;
for ind = 1:9
    % cycling from 1 to n-1 rather than 0 to n
    C = nchoosek(1:10,ind);
    P = ones(size(C,1),10);
    for indP = 1:size(C,1)
        P(indP,C(indP,:))=-1;
    end
    Permutations(count1 + (0:size(P,1)-1),:) = P;
    count1  = count1 + size(P,1);
end
Permutations(end+1,:) = repmat(-1,1,10); 
Permutations(end+1,:) = repmat(1,1,10);

%%  ==================== HG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,1,:,:)); % HG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(1,:) = p;
tval(1,:) = t_data;

resp_HG_unexp = datamat11;

%%  ==================== PP
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,2,:,:)); % PP

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(2,:) = p;
tval(2,:) = t_data;

resp_PP_unexp = datamat11;

%%  ==================== PT
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,3,:,:)); % PT

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat;  p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end
p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(3,:) = p;
tval(3,:) = t_data;

resp_PT_unexp = datamat11;

%%  ==================== aSTG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,4,:,:)); % aSTG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end
p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(4,:) = p;
tval(4,:) = t_data;

resp_aSTG_unexp = datamat11;

%%  ==================== pSTG
datamat1 = []; datamat11 = []; 
datamat1 = squeeze(datamat(:,:,5,:,:)); % pSTG

datamat11 = squeeze(nanmean(nanmean(datamat1,3),4)); % avg over conditions and hemi

[h,p,ci,stats] = ttest(datamat11); t_data = stats.tstat; p;

ts = [];
for i = 1:length(Permutations)
    [~,~,~,stats] = ttest(datamat11.*repmat(Permutations(i,:)',1,3));
    ts(i,:) = stats.tstat;
end

p = (sum(abs(ts)>= abs(t_data)))./length(Permutations); 

pvaluesCorr(5,:) = p;
tval(5,:) = t_data;

resp_pSTG_unexp = datamat11;

%% FDR correction

pval_vec = [pvaluesCorr(1,:), pvaluesCorr(2,:), pvaluesCorr(3,:), pvaluesCorr(4,:), pvaluesCorr(5,:)];

[FDR] = mafdr(pval_vec, 'BHFDR', 'true')       % Take all 15 numbers and correct.

pvaluesCorr = [FDR(1:3); FDR(4:6); FDR(7:9); FDR(10:12); FDR(13:15)];
pvaluesCorr = table(pvaluesCorr);
tvalTab = table(tval);
writetable(pvaluesCorr, ['pvalCorr', depths, '_permutations_FDRCorr_Unpred.txt'], 'Delimiter', '\t')
writetable(tvalTab, ['tvalTab', depths, '_permutations_FDRCorr_Unpred.txt'], 'Delimiter', '\t')


%% Save output file

save(['output_t2_model_Unpred', vox_tag, depths, '.mat'], 'resp_HG_unexp', 'resp_PP_unexp', 'resp_PT_unexp', ...
    'resp_aSTG_unexp', 'resp_pSTG_unexp');
