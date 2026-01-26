function [Tl_n, Th_n] = adaptiveCannyThresholdAuto(I)
% adaptiveCannyThresholdAuto
% Fully automatic Canny threshold selection (normalized)
% Bayesian + cross-entropy + posterior curvature
%
% Output:
%   Tl_n, Th_n : normalized low / high thresholds in [0,1]

%% Step 0: grayscale & double
if size(I,3) == 3
    I = rgb2gray(I);
end
I = double(I);

%% Step 1: gradient magnitude
[Gx, Gy] = imgradientxy(I);
G = sqrt(Gx.^2 + Gy.^2);
g = G(:);

%% Step 2: histogram (PDF)
numBins = 256;
[counts, edges] = histcounts(g, numBins);
binCenters = edges(1:end-1) + diff(edges)/2;
pdfG = counts / sum(counts);

%% Step 3: automatic prior_E (gradient tail)
medG = median(g);
madG = median(abs(g - medG));     % robust MAD
g0 = medG + 2 * madG;

prior_E = sum(g > g0) / numel(g);
prior_E = min(max(prior_E, 0.01), 0.3);

%% Step 4: non-parametric likelihoods
bgIdx   = binCenters <= g0;
edgeIdx = binCenters > g0;

pG_B = zeros(size(binCenters));
pG_E = zeros(size(binCenters));

pG_B(bgIdx)   = pdfG(bgIdx);
pG_E(edgeIdx) = pdfG(edgeIdx);

pG_B = pG_B / (sum(pG_B) + eps);
pG_E = pG_E / (sum(pG_E) + eps);

%% Step 5: Bayesian posterior
P_E = (pG_E * prior_E) ./ ...
      (pG_E * prior_E + pG_B * (1 - prior_E) + eps);

%% Step 6: High threshold (cross-entropy)
L = zeros(numBins,1);
for k = 1:numBins
    L(k) = -sum(log(1 - P_E(1:k) + eps)) ...
           -sum(log(P_E(k+1:end) + eps));
end
[~, idxTh] = min(L);
Th = binCenters(idxTh);

%% Step 7: Low threshold (posterior curvature / elbow)
P_s = smooth(P_E, 7);
d2P = abs(diff(P_s,2));
[~, idxTl] = max(d2P);
Tl = binCenters(idxTl);

%% Safety: ensure Tl < Th
if Tl >= Th
    Tl = 0.5 * Th;
end

%% Step 8: normalize to [0,1] for Canny
Gmax = max(g) + eps;
Tl_n = Tl / Gmax;
Th_n = Th / Gmax;

end
