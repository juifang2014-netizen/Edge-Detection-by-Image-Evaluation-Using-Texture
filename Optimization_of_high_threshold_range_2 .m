clear all
clc
close all
warning off

filepath='C:\Users\User\Desktop\新增資料夾';
cd(filepath);
codepath='C:\Users\User\Desktop\REVISE 1\程式\function';
files = dir(fullfile(filepath,'*.tif'));
if isempty(files)==1
    files = dir(fullfile(filepath,'*.png'));
end
if isempty(files)==1
    files = dir(fullfile(filepath,'*.tiff'));
end
lengthFiles = length(files);
% pause(180);

C_sigma=sqrt(1);
ThresholdRatio=0.1;
C_threshold_number=2^8;
C_threshold_number_threshold_0=0.95;
outlier_selection=1; % 去除 outluers: 1，不去除 outliers: 2
C_threshold_number_threshold=1-C_threshold_number_threshold_0;
d_ratio = 0.02;  % 容忍距離比例 (可選, 預設 0.02)
mode='all';  % 'all' (計算所有匹配平均距離, 預設) and 'TP'  (只計算 TP 匹配平均距離)
size_ratio=0.1;

size_order_threshold=1;

for i=1:lengthFiles;
    dlmwrite(['000_SizeRatio_' num2str(size_ratio) '.txt'], size_ratio);
    filesname_1_0 = dir(fullfile(filepath,'*.tif'));
    if isempty(filesname_1_0)==1
        filesname_1_0 = dir(fullfile(filepath,'*.png'));
    end
    if isempty(filesname_1_0)==1
        filesname_1_0 = dir(fullfile(filepath,'*.tiff'));
    end
    cell_str_1 = strsplit(filesname_1_0(i).name,'.');  % 分解字串，以 . 為界線
    filesname_1_1=cell_str_1{1};
    I_0 = imread(filesname_1_0(i).name);
    I=I_0;
    image_size=size(I);
    I_dimension=numel(image_size);
    if I_dimension==3
        if image_size(3)==4;
            I=I(:,:,1:3);
            I=rgb2gray(I);
        elseif image_size(3)==3
            I=rgb2gray(I);
        else
        end
    end
    I=im2double(I);
    
    %  Optimization of threshold dynamic range
    GaussianDieOff = .0001;
    pw = 1:30; % possible widths
    ssq = C_sigma^2;

    % 高斯平滑
    width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
    t = (-width:width);
    gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter
    [x,y]=meshgrid(-width:width,-width:width);
    dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
    aSmooth=imfilter(I,gau,'conv','replicate');   % run the filter across rows
    aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then across columns
    ix = imfilter(aSmooth, dgau2D, 'conv','replicate');
    iy = imfilter(aSmooth, dgau2D', 'conv','replicate');
    mag = sqrt((ix.*ix) + (iy.*iy));
    magmax = max(mag(:));
    if magmax > 0
        mag = mag / magmax;
    end
    
    dynamic_class_0=class(I_0);
    if strcmp(dynamic_class_0,'uint16')==1
        dynamic_class_0=2^16;
    elseif strcmp(dynamic_class_0,'uint8')==1
        dynamic_class_0=2^8;
    else
    end
    
    counts_1=imhist(mag, dynamic_class_0);    
    dynamic_class=2^16;
    C_high_threshold_0=1/dynamic_class:1/dynamic_class:1;
    counts_1=imresize(counts_1,dynamic_class/dynamic_class_0,'bilinear');
    counts_1=counts_1(:,1);
    
    % Method 2
    counts_2=counts_1;
    for LL=2:size(counts_1,1)-1
        counts_2(LL)=mean([counts_2(LL+1) counts_2(LL-1)]);
    end
    n_1_counts=counts_2/sum(counts_2);
    hs_counts=(-1)*n_1_counts.*(log2(n_1_counts));
    hs_counts=hs_counts(~isnan(hs_counts));  % 每個長條 in 直方圖
    hs_sum_counts(i)=sum(hs_counts);
    if min(find(counts_1==max(counts_1)))==1
        if hs_sum_counts(i)<11
            threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)/4);
        else
            threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)*2);
        end
    else
        threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)/2);
    end
%     threshold_high_0=threshold_range_0;
%     threshold_high_0(threshold_range_0<max(find(counts_2==max(counts_2))))=[];
%     threshold_high_0=min(threshold_high_0);    
%     TF_1_threshold=isempty(threshold_high_0);
%     if TF_1_threshold==1
%         threshold_high_0=size(counts_1,1);
%     else
%     end
    threshold_high_0=size(counts_1,1)*size_ratio;
    threshold_high_0=round(threshold_high_0);
    C_high_threshold_high=C_high_threshold_0(threshold_high_0);
    threshold_low_0=threshold_range_0;
    threshold_low_0(threshold_range_0>=max(find(counts_2==max(counts_2))))=[];
    threshold_low_0=min(threshold_low_0);
    TF_2_threshold=isempty(threshold_low_0);
    if  TF_2_threshold==1
        threshold_low_0=1;
    else
    end
    C_high_threshold_low=C_high_threshold_0(threshold_low_0);   
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(C_high_threshold_0,counts_1)
    hold on
    ylim([min(counts_1(:))-0.05*max(counts_1(:)) max(counts_1(:))+0.05*max(counts_1(:))])
    xlim([-0.025 1.025])
    scatter([C_high_threshold_high C_high_threshold_low],[counts_1(threshold_high_0) counts_1(threshold_low_0)],'r','filled')
    saveas(gcf,[filesname_1_1 '_high_threshold_of_dynamic_range.jpg']);  % 圖片存檔
    dlmwrite([filesname_1_1 '_intensity_histogram.txt'],counts_1);
    
    counts_1_first_half=counts_1(1:threshold_high_0);
    counts_1_second_half=counts_1(threshold_high_0:end);
    
    counts_1_ratio_first=sum(counts_1_first_half)/sum(counts_1);
    counts_1_ratio_first_1(i)=counts_1_ratio_first;
    close all
end

dlmwrite('Counts_ratio_first.txt',counts_1_ratio_first_1);

