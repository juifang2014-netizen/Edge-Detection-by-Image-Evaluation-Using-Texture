function result = computePRROC(labels, scores)
% PR + ROC metrics + diagnostic indicators (Best threshold, Sens, Spec, PPV, NPV, Accuracy)
    %% ---------- PR curve ----------
    [recall, precision, thPR] = perfcurve(labels, scores,1,'xCrit','reca','yCrit','prec');
    F1 = 2*(recall.*precision)./(recall+precision);
    [bestF1, idx] = max(F1);
    AUCPR = computeAUCPR(recall,precision);

    result.PR.recall    = recall;
    result.PR.precision = precision;
    result.PR.threshold = thPR;
    result.PR.F1        = F1;
    result.PR.bestF1    = bestF1;
    result.PR.bestRecall = recall(idx);
    result.PR.bestPrecision = precision(idx);
    result.PR.bestThreshold = thPR(idx);
    result.PR.AUCPR     = AUCPR;

    %% ---------- ROC curve ----------
    [FPR, TPR, thROC, ROCAUC] = perfcurve(labels, scores,1);

    Sens_all = TPR;              % Sensitivity at all thresholds
    Spec_all = 1-FPR;            % Specificity
    Acc_all  = (TPR*sum(labels) + Spec_all*sum(1-labels)) / length(labels);

    PPV_all  = (TPR*sum(labels)) ./ (TPR*sum(labels) + (1-Spec_all)*sum(1-labels));
    NPV_all  = (Spec_all*sum(1-labels)) ./ (Spec_all*sum(1-labels) + (1-TPR)*sum(labels));

    % Handle division by zero
    PPV_all(~isfinite(PPV_all)) = 0;
    NPV_all(~isfinite(NPV_all)) = 0;
    Acc_all(~isfinite(Acc_all))  = 0;

    % Best threshold: Youden index
    [~, idxBest] = max(Sens_all + Spec_all -1);

    result.ROC.FPR        = FPR;
    result.ROC.TPR        = TPR;
    result.ROC.threshold  = thROC;
    result.ROC.Specificity = Spec_all(idxBest);
    result.ROC.Sensitivity = Sens_all(idxBest);
    result.ROC.PPV         = PPV_all(idxBest);
    result.ROC.NPV         = NPV_all(idxBest);
    result.ROC.Accuracy    = Acc_all(idxBest);
    result.ROC.AUC         = ROCAUC;
    result.ROC.BestThreshold = thROC(idxBest);
end