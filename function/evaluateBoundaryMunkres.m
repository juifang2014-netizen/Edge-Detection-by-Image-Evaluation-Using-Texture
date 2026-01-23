function [Precision, Recall, F1, TP, FP, FN, assignment, meanError, totalCost] = evaluateBoundaryMunkres(GT, DE, imgSize, d_ratio, mode)
% evaluateBoundaryMunkres 使用 Munkres Algorithm 評估邊界偵測
%
% 輸入:
%   GT      - M x 2 Ground Truth 邊界點 (row, col)
%   DE      - N x 2 Detected Edges 邊界點 (row, col)
%   imgSize - [H W] 影像大小 (可選)
%   d_ratio - 容忍距離比例 (可選, 預設 0.02)
%   mode    - 'all' (計算所有匹配平均距離, 預設)
%             'TP'  (只計算 TP 匹配平均距離)
%
% 輸出:
%   Precision, Recall, F1 - 評估指標
%   TP, FP, FN            - 計算用的統計量
%   assignment            - Munkres 匹配向量
%   meanError             - 平均匹配距離
%   totalCost             - 所有匹配距離總和

if nargin < 3
    H = max([GT(:,1); DE(:,1)]);
    W = max([GT(:,2); DE(:,2)]);
else
    H = imgSize(1);
    W = imgSize(2);
end

if nargin < 4
    d_ratio = 0.02;  % 預設 2%
end

if nargin < 5
    mode = 'all';    % 預設計算所有匹配平均距離
end

%% 1. 建立成本矩陣
N = size(DE,1);
M = size(GT,1);
costMat = zeros(N, M);

for i = 1:N
    for j = 1:M
        costMat(i,j) = sqrt((DE(i,1)-GT(j,1))^2 + (DE(i,2)-GT(j,2))^2);
    end
end

%% 2. 設定容忍距離 d_max
d_max = round(d_ratio * max(H,W));

%% 3. 用 Munkres 做匹配
[assignment, ~] = munkres(costMat);

%% 4. 計算 TP, FP, FN, matchedDistances
TP = 0;
allDistances = []; % 所有匹配距離
TPDistances = [];  % TP 匹配距離

for i = 1:length(assignment)
    j = assignment(i);
    if j > 0
        dist = costMat(i,j);
        allDistances(end+1) = dist; %#ok<AGROW>
        if dist <= d_max
            TP = TP + 1;
            TPDistances(end+1) = dist; %#ok<AGROW>
        end
    end
end

FP = N - TP;
FN = M - TP;

%% 5. 計算 Precision, Recall, F1
Precision = TP / (TP + FP);
Recall = TP / (TP + FN);
F1 = 2 * Precision * Recall / (Precision + Recall);

%% 6. 計算 meanError 和 totalCost
totalCost = sum(allDistances);  % 所有匹配總成本

switch mode
    case 'all'
        if ~isempty(allDistances)
            meanError = mean(allDistances);
        else
            meanError = NaN;
        end
    case 'TP'
        if ~isempty(TPDistances)
            meanError = mean(TPDistances);
        else
            meanError = NaN;
        end
    otherwise
        error('mode must be ''all'' or ''TP''');
end

end
