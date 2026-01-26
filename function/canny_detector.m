function BW = canny_detector(I, sigma, threshold)
% I         : input image (gray or RGB)
% sigma     : Gaussian sigma
% threshold : [low high] in [0,1] or [] (auto, MATLAB-like)
%
% BW        : binary edge map

%% 0. Preprocess
if ndims(I) == 3
    I = rgb2gray(I);
end
I = single(I);        % 影像轉 single
I = I ./ max(I(:));   % normalize to [0,1]

%% 1. Derivative of Gaussian (DoG)
k = ceil(3*sigma);
[x,y] = meshgrid(-k:k, -k:k);

G  = exp(-(x.^2 + y.^2)/(2*sigma^2));
G  = G / sum(G(:));

Gx = -x .* G / sigma^2;
Gy = -y .* G / sigma^2;

% kernel 仍用 double
Ix = imfilter(I, double(Gx), 'replicate');
Iy = imfilter(I, double(Gy), 'replicate');

mag   = hypot(Ix, Iy);
theta = atan2(Iy, Ix) * 180 / pi;
theta(theta < 0) = theta(theta < 0) + 180;

%% 2. Vectorized Non-Maximum Suppression
nms = zeros(size(mag), 'single');

d0   = (theta < 22.5) | (theta >= 157.5);
d45  = (theta >= 22.5) & (theta < 67.5);
d90  = (theta >= 67.5) & (theta < 112.5);
d135 = (theta >= 112.5) & (theta < 157.5);

magE  = circshift(mag, [0 -1]);
magW  = circshift(mag, [0  1]);
magN  = circshift(mag, [-1 0]);
magS  = circshift(mag, [1  0]);
magNE = circshift(mag, [-1 -1]);
magNW = circshift(mag, [-1  1]);
magSE = circshift(mag, [1  -1]);
magSW = circshift(mag, [1   1]);

nms(d0)   = mag(d0)   .* (mag(d0)   >= magE(d0)  & mag(d0) >= magW(d0));
nms(d45)  = mag(d45)  .* (mag(d45)  >= magNE(d45)& mag(d45)>= magSW(d45));
nms(d90)  = mag(d90)  .* (mag(d90)  >= magN(d90) & mag(d90)>= magS(d90));
nms(d135) = mag(d135) .* (mag(d135) >= magNW(d135)& mag(d135)>= magSE(d135));

%% 3. Double Threshold
if nargin < 3 || isempty(threshold)
    % MATLAB-like automatic thresholds
    nz = nms(nms > 0);      
    nz = sort(nz(:));
    p = 0.7;               
    idx = max(1, round(p * numel(nz)));
    highT = nz(idx);
    lowT  = 0.4 * highT;
else
    maxN = max(nms(:));
    lowT  = threshold(1) * maxN;
    highT = threshold(2) * maxN;
end

strong = nms >= highT;
weak   = (nms >= lowT) & (nms < highT);

%% 4. Fast Hysteresis using connected components
BW = false(size(nms));
CC = bwconncomp(strong | weak, 8);

for k = 1:CC.NumObjects
    idx = CC.PixelIdxList{k};
    if any(strong(idx))
        BW(idx) = true;
    end
end
end