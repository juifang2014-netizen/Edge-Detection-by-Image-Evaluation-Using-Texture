clear all
clc
close all
warning off

%% ================================
% PRC & ROC FULL ANALYSIS (FINAL CLEAN)
%% ================================

%% -------------------------
% File settings
%% -------------------------
filepath = 'C:\Users\User\Desktop';
codepath = 'C:\Users\User\Desktop\程式';

features_name_label = 21;
features_name = { ...
    '1_entropy','2_skewness','3_kurtosis','4_autocorrelation', ...
    '5_cluster_prominence','6_cluster_shade','7_dissimilarity', ...
    '8_maximum_probability','9_angular_second_moment','10_contrast', ...
    '11_correlation','12_variance','13_inverse_difference_moment', ...
    '14_sum_average','15_sum_variance','16_sum_entropy', ...
    '17_entropy_haralick','18_difference_variance', ...
    '19_difference_entropy', ...
    '20_information_measure_of_correlation_1', ...
    '21_information_measure_of_correlation_2'};

inputFile   = fullfile(filepath,'1.xlsx');
summaryFile = fullfile(filepath,['PRC_ROC_summary_' features_name{features_name_label} '.xlsx']);
curveFile   = fullfile(filepath,['PRC_ROC_curve_points_' features_name{features_name_label} '.xlsx']);

%% -------------------------
% Load data
%% -------------------------
col = xlsread(inputFile);
posScore = col(:,1); posScore = posScore(~isnan(posScore));
negScore = col(:,2); negScore = negScore(~isnan(negScore));

scores = [posScore; negScore];
labels = [ones(length(posScore),1); zeros(length(negScore),1)];

%% -------------------------
% Compute PRC & ROC metrics
%% -------------------------
cd(codepath)
metrics = computePRROC(labels, scores);

%% -------------------------
% Bootstrap 95% CI
%% -------------------------
nBoot = 1000;
N = length(labels);
rng(0);

AUCPR_boot = nan(nBoot,1);
AUROC_boot = nan(nBoot,1);
Prec_boot  = nan(nBoot,1);
Rec_boot   = nan(nBoot,1);
F1_boot    = nan(nBoot,1);
Sens_boot  = nan(nBoot,1);
Spec_boot  = nan(nBoot,1);
Acc_boot   = nan(nBoot,1);

for b = 1:nBoot
    idx = randsample(N,N,true);
    yb = labels(idx);
    sb = scores(idx);

    if numel(unique(yb)) < 2
        continue
    end

    mb = computePRROCdiagnostic(yb,sb);

    AUCPR_boot(b) = mb.PR.AUCPR;
    AUROC_boot(b) = mb.ROC.AUC;
    Prec_boot(b)  = mb.PR.bestPrecision;
    Rec_boot(b)   = mb.PR.bestRecall;
    F1_boot(b)    = mb.PR.bestF1;
    Sens_boot(b)  = mb.ROC.Sensitivity;
    Spec_boot(b)  = mb.ROC.Specificity;
    Acc_boot(b)   = mb.ROC.Accuracy;
end

CI = @(x) prctile(x(~isnan(x)),[2.5 97.5]);

%% -------------------------
% Format output (4 decimals)
%% -------------------------
fmt = '%.4f (%.4f - %.4f)';

AUPRC_str = sprintf(fmt,metrics.PR.AUCPR,CI(AUCPR_boot));
AUROC_str = sprintf(fmt,metrics.ROC.AUC,CI(AUROC_boot));
Prec_str  = sprintf(fmt,metrics.PR.bestPrecision,CI(Prec_boot));
Rec_str   = sprintf(fmt,metrics.PR.bestRecall,CI(Rec_boot));
F1_str    = sprintf(fmt,metrics.PR.bestF1,CI(F1_boot));
Sens_str  = sprintf(fmt,metrics.ROC.Sensitivity,CI(Sens_boot));
Spec_str  = sprintf(fmt,metrics.ROC.Specificity,CI(Spec_boot));
Acc_str   = sprintf(fmt,metrics.ROC.Accuracy,CI(Acc_boot));

Th_PR_str  = sprintf('%.4f',metrics.PR.bestThreshold);
Th_ROC_str = sprintf('%.4f',metrics.ROC.BestThreshold);

%% -------------------------
% Summary table
%% -------------------------
SummaryHeader = { ...
    'Type','Threshold','AUPRC_or_AUROC', ...
    'Precision_or_Sensitivity','Recall_or_Specificity','F1_or_Accuracy'};

SummaryCell = {
    'PRC', Th_PR_str,  AUPRC_str, Prec_str, Rec_str, F1_str;
    'ROC', Th_ROC_str, AUROC_str, Sens_str, Spec_str, Acc_str
    };

xlswrite(summaryFile,[SummaryHeader; SummaryCell]);

%% -------------------------
% Auto-fit summary table columns
%% -------------------------
Excel = actxserver('Excel.Application');
Excel.Visible = false;

Workbook = Excel.Workbooks.Open(summaryFile);
Sheet = Workbook.Sheets.Item(1);
Sheet.Columns.AutoFit;

Workbook.Save;
Workbook.Close(false);
Excel.Quit;
delete(Excel);

%% =========================
% SAVE PRC / ROC CURVE POINTS (ONLY 2 SHEETS)
%% =========================

% 刪除舊檔案
if exist(curveFile,'file')
    delete(curveFile);
end

% ---------- PRC ----------
minLen = min([ ...
    numel(metrics.PR.threshold), ...
    numel(metrics.PR.recall), ...
    numel(metrics.PR.precision), ...
    numel(metrics.PR.F1)]);

PRC_cell = [ ...
    {'Threshold','Recall','Precision','F1'}; ...
    num2cell([ ...
        metrics.PR.threshold(1:minLen), ...
        metrics.PR.recall(1:minLen), ...
        metrics.PR.precision(1:minLen), ...
        metrics.PR.F1(1:minLen)])];

xlswrite(curveFile, PRC_cell, 'PRC_curve');

% ---------- ROC ----------
[Xroc,Yroc,Troc] = perfcurve(labels,scores,1);

ROC_cell = [ ...
    {'Threshold','FPR','TPR'}; ...
    num2cell([Troc, Xroc, Yroc])];

xlswrite(curveFile, ROC_cell, 'ROC_curve');

%% -------------------------
% Plot PRC
%% -------------------------
figure;
plot(metrics.PR.recall,metrics.PR.precision,'b','LineWidth',2); hold on
scatter(metrics.PR.bestRecall,metrics.PR.bestPrecision,80,'r','filled')
xlabel('Recall'); ylabel('Precision')
title(sprintf('PRC (AUPRC = %.4f)',metrics.PR.AUCPR))
grid on
print(gcf,fullfile(filepath,['PRC_' features_name{features_name_label} '.jpg']),'-djpeg','-r300')

%% -------------------------
% Plot ROC
%% -------------------------
figure;
plot(Xroc,Yroc,'k','LineWidth',2); hold on
plot([0 1],[0 1],'k--')
[~,idx] = min(abs(Troc-metrics.ROC.BestThreshold));
scatter(Xroc(idx),Yroc(idx),80,'r','filled')
xlabel('False Positive Rate'); ylabel('True Positive Rate')
title(sprintf('ROC (AUROC = %.4f)',metrics.ROC.AUC))
grid on
print(gcf,fullfile(filepath,['ROC_' features_name{features_name_label} '.jpg']),'-djpeg','-r300')

cd(filepath)
close all
