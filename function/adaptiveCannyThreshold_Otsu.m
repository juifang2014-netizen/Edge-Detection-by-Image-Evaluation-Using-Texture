function [Th_norm, Tl_norm] = adaptiveCannyThresholdPlot(I, lowHighRatio)
% adaptiveCannyThresholdPlot 計算自適應 Canny 高低閾值並繪圖
%   [Th_norm, Tl_norm] = adaptiveCannyThresholdPlot(I, lowHighRatio)
%   I: 輸入影像（灰階或彩色）
%   lowHighRatio: 可選，低閾值 / 高閾值比例（若不給則自適應）
%
%   Th_norm: 高閾值（0~1，可直接給 edge 使用）
%   Tl_norm: 低閾值（0~1，可直接給 edge 使用）

    if nargin < 2
        lowHighRatio = []; % 空值表示自適應
    end

    % 如果是彩色影像，轉灰階
    if size(I,3) == 3
        I = rgb2gray(I);
    end
    I = im2double(I);

    %% 1. 計算梯度 (Sobel)
    [Gx, Gy] = imgradientxy(I, 'sobel');
    [Gmag, Gdir] = imgradient(Gx, Gy);

    %% 2. 非極大值抑制
    G_nms = nonMaxSuppression(Gmag, Gdir);

    %% 3. 計算高閾值
    G_vec = G_nms(:);
    G_vec = G_vec(G_vec > 0);  % 移除背景
    if isempty(G_vec)
        Th_norm = 0;
        Tl_norm = 0;
        warning('No gradient detected.');
        return;
    end
    Gmax = max(G_vec);           % 最大梯度
    G_norm = G_vec / Gmax;       % 正規化到 [0,1]

    % 使用 OTSU 計算高閾值
    level = graythresh(G_norm);  
    Th_norm = level;             % 高閾值 0~1

    %% 4. 計算低閾值
    if isempty(lowHighRatio)
        % 自適應低閾值比率，完全不限制範圍
        muG = mean(G_vec);
        sigmaG = std(G_vec);
        ratio_adaptive = muG / (muG + sigmaG); % 梯度集中 -> ratio 高，分散 -> ratio 低
        Tl_norm = ratio_adaptive * Th_norm;
    else
        Tl_norm = lowHighRatio * Th_norm;
    end

    %% 5. 繪製梯度直方圖 + 高低閾值
%     figure;
%     histogram(G_vec, 100, 'FaceColor', [0.2 0.6 0.8]); hold on;
%     ylims = ylim;
%     plot([Th_norm*Gmax Th_norm*Gmax], ylims, 'r', 'LineWidth', 2, 'DisplayName', 'High threshold (Th)');
%     plot([Tl_norm*Gmax Tl_norm*Gmax], ylims, 'g', 'LineWidth', 2, 'DisplayName', 'Low threshold (Tl)');
%     xlabel('Gradient magnitude');
%     ylabel('Pixel count');
%     title('Gradient histogram with adaptive thresholds');
%     legend('show');
%     grid on;

end
