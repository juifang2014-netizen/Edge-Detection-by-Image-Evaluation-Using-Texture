function AUCPR = computeAUCPR(recall, precision)
    recall = recall(:); precision = precision(:);
    valid = isfinite(recall) & isfinite(precision);
    recall = recall(valid); precision = precision(valid);
    if numel(recall)<2, AUCPR=NaN; return; end
    [recall, idx] = sort(recall); precision = precision(idx);
    [recall, ia] = unique(recall,'stable'); precision = precision(ia);
    AUCPR = trapz(recall,precision);
end