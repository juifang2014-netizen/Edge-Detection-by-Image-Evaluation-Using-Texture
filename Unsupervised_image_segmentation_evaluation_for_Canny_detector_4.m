clear all
clc
close all
warning off

filepath='C:\Users\User\Desktop\新增資料夾';
codepath='C:\Users\User\Desktop\REVISE 1\程式\function';
% filename=length(dir('*.tif'));  % 資料夾內檔案數
files = dir(fullfile(filepath,'*.tif'));
if isempty(files)==1
    files = dir(fullfile(filepath,'*.png'));
end
if isempty(files)==1
    files = dir(fullfile(filepath,'*.tiff'));
end
lengthFiles = length(files);
% pause(120);

C_sigma=sqrt(1);
% C_sigma=5;
% C_low_threshold=0.02;
% C_high_threshold=0.03;
% C_high_threshold=2.5*C_low_threshold;
% Gaussian_sigma=30;
% PercentOfPixelsNotEdges=0.7;
ThresholdRatio=0.01;
C_threshold_number=2^8;
C_threshold_number_threshold_0=0.95;
outlier_selection=1; % 去除 outluers: 1，不去除 outliers: 2
% threshold_range_selection=2; % 0.5: 1，2: 2
C_threshold_number_threshold=1-C_threshold_number_threshold_0;
d_ratio = 0.02;  % 容忍距離比例 (可選, 預設 0.02)
mode='all';  % 'all' (計算所有匹配平均距離, 預設) and 'TP'  (只計算 TP 匹配平均距離)

size_order_threshold=1;
% C_threshold_selection=1;
% if C_threshold_selection==1
%     C_low_threshold_1=0.001:0.00001:0.01;
% else
%     C_low_threshold_1=0.01:0.001:0.1;
%     C_low_threshold_1=0.001:0.0001:0.02;
% end

for i=1:lengthFiles;
    % 每張圖的 threshold range 不同，需要清除相關參數
    clear A AA threshold_range_0 C_high_threshold C_high_threshold_1 distance_1 ROI_parameters Background_parameters Dice_1 Jaccard_1 Dice_Jaccard_index Accuracy Sensitivity Fmeasure Precision MCC Dice Jaccard Specitivity cf_points_sigmoid_4_parameters fx_sigmoid_4_parameters_1 fxx_sigmoid_4_parameters_1 distance_fitting_r2_sigmoid_4_parameters location_min_fitting_slope_sigmoid_4_parameters location_mirror_min_fitting_slope_sigmoid_4_parameters cf_points_poly fx_poly_1 fxx_poly_1 distance_fitting_r2_poly location_appropriate_intersection_distance_fitting location_optimal
    tic;
    kk=1;
    cd(filepath);
    dlmwrite(['000_ThresholdRatio_' num2str(ThresholdRatio) '.txt'], ThresholdRatio);
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
%     I_0=imcomplement(I_0);
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
%     I=histeq(I);
    
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
%     counts_1=counts_1(:,1);
%     counts_1=resample(counts_1,dynamic_class/dynamic_class_0,1);
%     counts_1=interp(counts_1,dynamic_class/dynamic_class_0);
    
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
% %     imhist(mag, dynamic_class)
%     plot(C_high_threshold_0,counts_1)
%     hold on
%     ylim([0 max(counts_1(:))+0.05*max(counts_1(:))])
%     xlim([-0.025 1])
    
    % Method 2
    counts_2=counts_1;
    for LL=2:size(counts_1,1)-1
        counts_2(LL)=mean([counts_2(LL+1) counts_2(LL-1)]);
    end
%     A=graycomatrix(mag);
%     AA=graycoprops(A,{'contrast'});
%     AAA(i)=AA.Contrast;
%     if AAA(i)<0.05
%         threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)/2);
%     else
%         threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)*2);
%     end
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
%     threshold_range_0=find(counts_2<=((size(I,1)*size(I,2))/dynamic_class_0)/2);
%     threshold_range_0=find(counts_2<=max(counts_1-3*std(counts_1)));
    threshold_high_0=threshold_range_0;
    threshold_high_0(threshold_range_0<max(find(counts_2==max(counts_2))))=[];
    threshold_high_0=min(threshold_high_0);    
%     if threshold_high_0>size(counts_1,1)
%         threshold_high_0=size(counts_1,1);
%     else
%     end
    TF_1_threshold=isempty(threshold_high_0);
    if TF_1_threshold==1
        threshold_high_0=size(counts_1,1);
    else
    end
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
%     imhist(mag, dynamic_class)
    plot(C_high_threshold_0,counts_1)
    hold on
    ylim([min(counts_1(:))-0.05*max(counts_1(:)) max(counts_1(:))+0.05*max(counts_1(:))])
    xlim([-0.025 1.025])
%     scatter(1/dynamic_class,counts_1(1),'c','filled');
%     scatter(threshold_high,counts_1(threshold_high_0),'r','filled');
    scatter([C_high_threshold_high C_high_threshold_low],[counts_1(threshold_high_0) counts_1(threshold_low_0)],'r','filled')
    saveas(gcf,[filesname_1_1 '_high_threshold_of_dynamic_range.jpg']);  % 圖片存檔
    dlmwrite([filesname_1_1 '_intensity_histogram.txt'],counts_1);
   
    % Optimization of threshold numbers
    % Method 2
    for ii=1:threshold_high_0
        m_mean(ii)=mean(counts_1(threshold_low_0:ii:threshold_high_0));
        m_medi(ii)=median(counts_1(threshold_low_0:ii:threshold_high_0));
        m_max(ii)=max(counts_1(threshold_low_0:ii:threshold_high_0));
%         m_max(ii)=mean(counts_1(threshold_low_0:ii:threshold_high_0));
        s_std(ii)=std(counts_1(threshold_low_0:ii:threshold_high_0));
%         s_std(ii)=rms(counts_1(threshold_low_0:ii:threshold_high_0));
        k_kurt(ii)=kurtosis(counts_1(threshold_low_0:ii:threshold_high_0));
        s_skew(ii)=skewness(counts_1(threshold_low_0:ii:threshold_high_0));
    end
    for ii_mean=1:threshold_high_0
        if m_mean(ii_mean)<m_mean(1)-C_threshold_number_threshold*m_mean(1) || m_mean(ii_mean)>m_mean(1)+C_threshold_number_threshold*m_mean(1);
            break
        end
    end
    for ii_medi=1:threshold_high_0
        if m_medi(ii_medi)<m_medi(1)-C_threshold_number_threshold*m_medi(1) || m_medi(ii_medi)>m_medi(1)+C_threshold_number_threshold*m_medi(1);
            break
        end
    end
    for ii_max=1:threshold_high_0
        if m_max(ii_max)<m_max(1)-C_threshold_number_threshold*m_max(1) || m_max(ii_max)>m_max(1)+C_threshold_number_threshold*m_max(1);
            break
        end
    end
    for ii_std=1:threshold_high_0
        if s_std(ii_std)<s_std(1)-C_threshold_number_threshold*s_std(1) || s_std(ii_std)>s_std(1)+C_threshold_number_threshold*s_std(1);
            break
        end
    end
    for ii_kurt=1:threshold_high_0
        if k_kurt(ii_kurt)<k_kurt(1)-C_threshold_number_threshold*k_kurt(1) || k_kurt(ii_kurt)>k_kurt(1)+C_threshold_number_threshold*k_kurt(1);
            break
        end
    end
    for ii_skew=1:threshold_high_0
        if s_skew(ii_skew)<s_skew(1)-C_threshold_number_threshold*s_skew(1) || s_skew(ii_skew)>s_skew(1)+C_threshold_number_threshold*s_skew(1);
            break
        end
    end

%     if dynamic_class_0==2^16;
%         C_threshold_number=max(ii_max,ii_std);
%     elseif dynamic_class_0==2^8
%         C_threshold_number=min(ii_max,ii_std);
%     else
%     end   
%     C_threshold_number=min(ii_max-1,ii_std-1);
% 
%     % Method 3
%     for ii=1:threshold_high_0
%         m_max(ii)=max(counts_1(threshold_low_0:ii:threshold_high_0));
%         s_std(ii)=std(counts_1(threshold_low_0:ii:threshold_high_0));
%     end
%     ii_max=findchangepts(m_max,'Statistic','linear','MaxNumChanges',2);
%     ii_max=min(ii_max);
%     ii_std=findchangepts(s_std,'Statistic','linear','MaxNumChanges',2);
%     ii_std=min(ii_std);

%     C_threshold_number=min([ii_mean ii_medi ii_max ii_std ii_kurt ii_skew]);
%     C_threshold_number=median([ii_mean ii_medi ii_max ii_std ii_kurt ii_skew]);
    C_threshold_number=prctile([ii_mean ii_medi ii_max ii_std ii_kurt ii_skew],50);
    dlmwrite([filesname_1_1 '_threshold_number.txt'],[C_threshold_number ii_mean ii_medi ii_max ii_std ii_kurt ii_skew]);
    dlmwrite([filesname_1_1 '_high_threshold_moving_mean_median_max_standard_deviation_kurttosis_skewness.txt'],[m_mean' m_medi' m_max' s_std' k_kurt' s_skew']);
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(m_mean)
    hold on
    scatter(ii_mean,m_mean(ii_mean),'r','filled');
    saveas(gcf,[filesname_1_1 '_mean_and_threshold_number.jpg']);  % 圖片存檔
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(m_medi)
    hold on
    scatter(ii_medi,m_medi(ii_medi),'r','filled');
    saveas(gcf,[filesname_1_1 '_median_and_threshold_number.jpg']);  % 圖片存檔
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(m_max)
    hold on
    scatter(ii_max,m_max(ii_max),'r','filled');
    saveas(gcf,[filesname_1_1 '_max_and_threshold_number.jpg']);  % 圖片存檔
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(s_std)
    hold on
    scatter(ii_std,s_std(ii_std),'r','filled');
    saveas(gcf,[filesname_1_1 '_standard_deviation_and_threshold_number.jpg']);  % 圖片存檔
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(k_kurt)
    hold on
    scatter(ii_kurt,k_kurt(ii_kurt),'r','filled');
    saveas(gcf,[filesname_1_1 '_kurtosis_and_threshold_number.jpg']);  % 圖片存檔
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    plot(s_skew)
    hold on
    scatter(ii_skew,s_skew(ii_skew),'r','filled');
    saveas(gcf,[filesname_1_1 '_skewness_and_threshold_number.jpg']);  % 圖片存檔
    C_high_threshold_1=C_high_threshold_low:1/dynamic_class*C_threshold_number:C_high_threshold_high;    
   
%     C_high_threshold_1=1/dynamic_class:(1/dynamic_class)*(dynamic_class/256):threshold_high;    
    
%     dlmwrite(['fitting_points_high_threshold_range_'  filesname_1_1 '_' num2str(C_sigma) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '_exp1.txt'],cf_exp1_points_0);
%     dlmwrite(['fitting_r2_high_threshold_range_'  filesname_1_1 '_' num2str(C_sigma) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '_exp1.txt'],G_exp1.rsquare);
%     dlmwrite([filesname_1_1 '_threshold_range.txt'],[ipt_high C_high_threshold_high ipt_low C_high_threshold_low]);
    dlmwrite([filesname_1_1 '_threshold_range.txt'],[threshold_high_0 C_high_threshold_high threshold_low_0 C_high_threshold_low]);
    dlmwrite(['C_high_threshold_'  filesname_1_1 '_' num2str(C_sigma^2) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '.txt'],C_high_threshold_1');
    
%     C_high_threshold_1=min(mag(:)):C_threshold_interval/dynamic_class:max(mag(:));
%     dynamic_class=max(max(I_0))-min(min(I_0))+1;
%     I_1=I;
% %     I_1(I_1>mean(I_1(:))+3*std(I_1(:)))=0;

%     threshold_dynamic_range_max=max(counts_1);
%     threshold_dynamic_range_max_location=find(counts_1==threshold_dynamic_range_max);
%     threshold_dynamic_range_max_location_1=1/dynamic_class*threshold_dynamic_range_max_location;
%     counts_1_zero_location=find(counts_1==0);
%     counts_2=counts_1;  
%         
%     for LL_zero=2:size(counts_1_zero_location,1)-1
%         counts_2(counts_1_zero_location(LL_zero))=mean([counts_1(counts_1_zero_location(LL_zero)+1) counts_1(counts_1_zero_location(LL_zero)-1)]);
%     end    
%     for LL_high=threshold_dynamic_range_max_location:dynamic_class
%         if counts_2(LL_high)==0
%             break
%         end
%     end
%      for LL_low=threshold_dynamic_range_max_location:-1:1
%         if counts_2(LL_low)==0
%             break
%         end
%      end
%        
%     threshold_high=1/dynamic_class*LL_high;
%     threshold_low=1/dynamic_class*LL_low;
%     C_high_threshold_1=threshold_low:C_threshold_interval/dynamic_class:threshold_high-1/dynamic_class;
%     C_threshold_interval=C_threshold_multi;
%     C_high_threshold_1=LL_low:C_threshold_interval:LL_high;
    
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%     imhist(mag,dynamic_class);
%     hold on
%     scatter(threshold_high,counts_1(LL_high),'r','filled')
%     scatter(threshold_low,counts_1(LL_low),'r','filled')
%     xlim([min(mag(:)) max(mag(:))])
%     ylim([0 max(counts_1)+0.05*max(counts_1)])
%     plot(counts_1);
%     hold on
%     scatter(LL_high,counts_1(LL_high),'r','filled')
%     scatter(LL_low,counts_1(LL_low),'r','filled')
%     saveas(gcf,[filesname_1_1 '_threshold_dynamic_range.jpg']);  % 圖片存檔
%     dlmwrite([filesname_1_1 '_thrshold_dynamic_range.txt'], [LL_high threshold_high  LL_low threshold_low]);

%     [BW_1,thres]=edge(I,'Canny');
%     cd(codepath);
%     [lowThresh, highThresh] = Canny_threshold(I,C_sigma,PercentOfPixelsNotEdges,ThresholdRatio);
%     cd(filepath);
%     C_high_threshold_1=(C_threshold_interval/PercentOfPixelsNotEdges)*highThresh:(C_threshold_interval/PercentOfPixelsNotEdges)*highThresh:(1/C_threshold_interval)*(C_threshold_interval/PercentOfPixelsNotEdges)*highThresh;

    bw_manual_segmentation=load(['RANGE_' filesname_1_1 '.txt']);  % 輸入手動分割的結果
    
%     for C_low_threshold=C_low_threshold_1;
%         C_high_threshold=0.0001+C_low_threshold; 

%     for C_high_threshold=C_high_threshold_1
%         C_low_threshold_1=C_high_threshold*ThresholdRatio;
%         C_low_threshold=C_low_threshold_1;

%     for PercentOfPixelsNotEdges=threshold_low:C_threshold_interval:threshold_high

    for C_high_threshold=C_high_threshold_1
        C_low_threshold=C_high_threshold*ThresholdRatio;
%         cd(codepath);
%         [C_low_threshold, C_high_threshold] = Canny_threshold(I,C_sigma,PercentOfPixelsNotEdges,ThresholdRatio, dynamic_class);
%         cd(filepath);
%         filesname_1_0 = dir(fullfile(filepath,'*.tif'));
%         cell_str_1 = strsplit(filesname_1_0(i).name,'.');  % 分解字串，以 . 為界線
%         filesname_1_1=cell_str_1{1};
        filesname_1=['_' num2str(C_sigma^2) '_' num2str(C_high_threshold) '_' num2str(C_low_threshold)];
%         I_0 = imread(filesname_1_0(i).name);
%         I=I_0;
%         I=im2double(I);
%         C_low_threshold_1(kk)=C_low_threshold;
%         C_high_threshold_1(kk)=C_high_threshold;        
        
%         BW = edge(I,'Canny',[C_low_threshold C_high_threshold], C_sigma);
        cd(codepath);
        BW = canny_detector(I, C_sigma, [C_low_threshold C_high_threshold]);
        cd(filepath);
        
        sel1=strel('square',3);
        sel2=strel('square',5);
        BW=imdilate(BW,sel2);
%         BW=imerode(BW,sel1);
        BW = imfill(BW, 'holes');

%         BW=imgaussfilt(uint8(BW),Gaussian_sigma);
% 
%         BW=imcomplement(BW);
        imLabel = bwlabel(BW);% ??通?域?行?
        stats = regionprops(imLabel,'Area');
        [b,index]=sort([stats.Area],'descend');
        if length(stats)<size_order_threshold  % 前幾大
            bw=imLabel;
        else
            bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
        end
        
        % Performance evaluation of edge edtection
        cd(codepath);
        [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
        cd(filepath);
        Dice_1(kk)=Dice;
        Jaccard_1(kk)=Jaccard;
%         Dice_Jaccard_index=[Dice_1;Jaccard_1]';
        
%         bw = imfill(bw, 'holes'); 
%         dlmwrite([filesname_1 '_ROI.txt'],bw);
        bw_edge=bwboundaries(bw);
        k=size(bw_edge,1);
        thisBoundary = bw_edge{k};
        thisBoundary_y=thisBoundary(:,2);
        thisBoundary_x=thisBoundary(:,1);
%         dlmwrite([filesname_1 '_thisBoundary_y.txt'],thisBoundary_y);
%         dlmwrite([filesname_1 '_thisBoundary_x.txt'],thisBoundary_x);
        figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
        I_1=I_0;
        if I_dimension==3
            if image_size(3)==4;
                I_1=I_1(:,:,1:3);
                if dynamic_class_0==2^16;
                    I_1(I_1>=2500)=2500;  % for brightfield image
                    imagesc(I_1)
                     colormap(gray)
                else
                    imshow(I_1)
                end
            elseif image_size(3)==3;
                if dynamic_class_0==2^16;
                    I_1(I_1>=2500)=2500;  % for brightfield image
                    imagesc(I_1)
                     colormap(gray)
                else
                    imshow(I_1)
                end
            end
        else
            if dynamic_class_0==2^16;
                I_1(I_1>=2500)=2500;  % for brightfield image
                imagesc(I_1)
                colormap(gray)
            else
                imshow(I_1)
            end
        end
        axis image
        axis off
        saveas(gcf,[filesname_1_1 '_1.jpg']);  % 圖片存檔
        hold on
        plot(thisBoundary_y, thisBoundary_x,'b:', 'LineWidth', 3)
        saveas(gcf,[filesname_1_1 '_' num2str(kk) filesname_1 '_with_ROI.jpg']);  % 圖片存檔
%         [xi,xj]=find(I_0>2500);
%         imagesc(I_1)
%         colormap(gray)
%         hold on
%         scatter(xj,xi,'r','filled')

        % Performance evaluation of edge edtection
        p=load(['ROI_' filesname_1_1 '.txt']);
        q=[thisBoundary_y thisBoundary_x];
        cd(codepath);
        p = resampleBoundary(p, 200);
        q = resampleBoundary(q, 200);        
        [Precision, Recall, F1, TP, FP, FN, assignment, meanError, totalCost] = evaluateBoundaryMunkres(p, q, size(I), d_ratio, mode);
        cd(filepath)
        Precision_1(kk)=Precision;
        meanError_1(kk)=meanError;
        totalCost_1(kk)=totalCost;
        
        % ROI parameter
%         % Shape parameter
%         bw_1_ROI = edge(bw);
%         % Fractal dimension
%         cd(codepath);
%         [ D_ROI ] = hausDim( bw_1_ROI );
%         D_1_ROI(kk)=D_ROI;
%         cd(filepath);
%         % Basic region parameters
%         stats_regionprops_ROI=regionprops(bw,'Area','Centroid','Perimeter','Eccentricity','Extent','Extrema','Orientation','MajorAxisLength','MinorAxisLength');
%         % Circularity
%         allAreas_ROI(kk) = [stats_regionprops_ROI.Area];
%         allPerimeters_ROI(kk) = [stats_regionprops_ROI.Perimeter];
%         circularity_ROI(kk) = (4 * pi * allAreas_ROI(kk))/(allPerimeters_ROI(kk)^ 2);
        
        % Texture parameter
        I_ROI=I.*logical(bw);
%         I_ROI_1=I_ROI;
%         I_ROI_1=I_ROI_1(:);
        ROI_intensity_0=I_ROI;
%         if mean(mean(ROI_intensity_0))==0
%             ROI_intensity_0(1:3,1:3)=min(I_ROI_1(find( I_ROI_1-min( I_ROI_1))));            
%         else
%         end
        ROI_intensity_0(ROI_intensity_0==0)=[];
        % Intensity
%         ROI_intensity=mean(mean(ROI_intensity_0));
        % 標準化(Z 分數)
        ROI_intensity_0=zscore(ROI_intensity_0);
        % Shannon entropy
        histnumber=32;
        figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
        n_ROI=histogram(ROI_intensity_0,histnumber);
%         saveas(gcf,[filesname_1 '_ROI_histogram_' num2str(histnumber) '.jpg']);  % 圖片存檔
        n_2_ROI=get(n_ROI,'Values');  % 讀取 histogram 相關數據時不能關圖
        n_1_ROI=n_2_ROI/size(ROI_intensity_0,2); 
        hs_ROI=(-1)*n_1_ROI.*(log2(n_1_ROI));
        hs_ROI=hs_ROI(~isnan(hs_ROI));  % 每個長條 in 直方圖
        hs_sum_ROI=sum(hs_ROI);   % 每個直方圖的範圍
        % Skwness
        s_ROI=skewness(ROI_intensity_0);
        % kurtosis
        k_texture_ROI=kurtosis(ROI_intensity_0);
        % GLCM
        glcms_ROI = graycomatrix(ROI_intensity_0);
%         stats_GLCM_ROI = graycoprops(glcms_ROI,{'contrast','Homogeneity','Correlation','Energy'});
%         Con_ROI=stats_GLCM_ROI.Contrast;
%         Homo_ROI=stats_GLCM_ROI.Homogeneity;
%         Corr_ROI=stats_GLCM_ROI.Correlation;
%         En_ROI=stats_GLCM_ROI.Energy;
        % Haralick features
        cd(codepath);
        [Haralick_ROI] = haralickTextureFeatures(glcms_ROI);
        [out_ROI] = GLCM_Features1(glcms_ROI);
        cd(filepath);
%         ROI_parameters_00=[ROI_intensity;hs_sum_ROI;s_ROI;k_texture_ROI;Con_ROI;Homo_ROI;Corr_ROI;En_ROI];
        ROI_parameters_00=[hs_sum_ROI;s_ROI;k_texture_ROI];
        ROI_parameters_0=[ROI_parameters_00;out_ROI.autoc;out_ROI.cprom;out_ROI.cshad;out_ROI.dissi;out_ROI.maxpr;Haralick_ROI(1:13)];
        ROI_parameters(:,kk)=ROI_parameters_0;
%         ROI_parameters(:,kk)=[ROI_intensity;hs_sum_ROI;s_ROI;k_texture_ROI];

        % Background parameter
%         % Shape parameter
%         bw_1_Background = edge(imcomplement(bw));
%         % Fractal dimension
%         cd(codepath);
%         [ D_Background ] = hausDim( bw_1_Background );
%         D_1_Background(kk)=D_Background;
%         cd(filepath);
%         % Basic region parameters
%         stats_regionprops_Background=regionprops(imcomplement(bw),'Area','Centroid','Perimeter','Eccentricity','Extent','Extrema','Orientation','MajorAxisLength','MinorAxisLength');
%         % Circularity
%         allAreas_Background(kk) = [stats_regionprops_Background.Area];
%         allPerimeters_Background(kk) = [stats_regionprops_Background.Perimeter];
%         circularity_Background(kk) = (4 * pi * allAreas_Background(kk))/(allPerimeters_Background(kk)^ 2);
        
        % Texture parameter
        I_background=I.*logical(imcomplement(bw));
        Background_intensity_0=I_background;
%         if mean(mean(Background_intensity_0))==0
%             Background_intensity_0(1:3,1:3)=min(I_ROI_1(find( I_ROI_1-min( I_ROI_1))));            
%         end
        Background_intensity_0(Background_intensity_0==0)=[];        
        % Intensity
%         Background_intensity=mean(mean(Background_intensity_0));
        % 標準化(Z 分數)
        Background_intensity_0=zscore(Background_intensity_0);        
        % Shannon entropy
        histnumber=32;
        figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
        n_background=histogram(Background_intensity_0,histnumber);
%         saveas(gcf,[filesname_1 '_background_histogram_' num2str(histnumber) '.jpg']);  % 圖片存檔
        n_2_background=get(n_background,'Values');  % 讀取 histogram 相關數據時不能關圖
        n_1_background=n_2_background/size(Background_intensity_0,2); 
        hs_background=(-1)*n_1_background.*(log2(n_1_background));
        hs_background=hs_background(~isnan(hs_background));  % 每個長條 in 直方圖
        hs_sum_background=sum(hs_background);   % 每個直方圖的範圍
        % Skwness
        s_background=skewness(Background_intensity_0);
        % kurtosis
        k_texture_background=kurtosis(Background_intensity_0);
        % GLCM
        glcms_Background = graycomatrix(Background_intensity_0);
%         stats_GLCM_Background = graycoprops(glcms_Background,{'contrast','Homogeneity','Correlation','Energy'});
%         Con_Background=stats_GLCM_Background.Contrast;
%         Homo_Background=stats_GLCM_Background.Homogeneity;
%         Corr_Background=stats_GLCM_Background.Correlation;
%         En_Background=stats_GLCM_Background.Energy;
        % Haralick features
        cd(codepath);
        [Haralick_Background] = haralickTextureFeatures(glcms_Background);
        [out_Background] = GLCM_Features1(glcms_Background);
        cd(filepath);
%         Background_parameters_00=[Background_intensity;hs_sum_background;s_background;k_texture_background;Con_Background;Homo_Background;Corr_Background;En_Background];
        Background_parameters_00=[hs_sum_background;s_background;k_texture_background];
        Background_parameters_0=[Background_parameters_00;out_Background.autoc;out_Background.cprom;out_Background.cshad;out_Background.dissi;out_Background.maxpr;Haralick_Background(1:13)];
        Background_parameters(:,kk)=Background_parameters_0;
%         Background_parameters(:,kk)=[Background_intensity;hs_sum_background;s_background;k_texture_background];

%         % Euclidean distance        
%         euclidean_distance_0=([ROI_parameters(:,kk) Background_parameters(:,kk)])';
%         euclidean_distance(kk)=pdist(euclidean_distance_0,'euclidean');        
       
        kk=kk+1;
        close all
%         estimation_time(i)=toc;          
    end
    
%     dlmwrite(['Estimation_time_' filesname_1_1 '.txt'],estimation_time);
    distance_1=sqrt((ROI_parameters'-Background_parameters').^2);
    %  Canberra distance 
    for dd=1:size(ROI_parameters,2)
        ppp=ROI_parameters(:,dd);
        qqq=Background_parameters(:,dd);
        distance_1(dd)=sum(abs(ppp - qqq)./(abs(ppp) + abs(qqq)));
    end 
           
    dlmwrite(['ROI_parameters_'  filesname_1_1 '.txt'],ROI_parameters');
    dlmwrite(['Background_parameters_'  filesname_1_1 '.txt'],Background_parameters');
    dlmwrite(['Distance_'  filesname_1_1 '.txt'],distance_1);
    dlmwrite(['Jaccard_index_'  filesname_1_1 '.txt'],Jaccard_1);
    dlmwrite(['Dice_index_'  filesname_1_1 '.txt'],Dice_1);
    dlmwrite(['Precision_'  filesname_1_1 '.txt'],Precision_1);
    dlmwrite(['meanError_'  filesname_1_1 '.txt'],meanError_1);
    dlmwrite(['totalCost_'  filesname_1_1 '.txt'],totalCost_1);
    Dice_index_max(i,:)=[max(Dice_1) min(find(Dice_1==max(Dice_1))) C_high_threshold_1(min(find(Dice_1==max(Dice_1))))]';
    Jaccard_index_max(i,:)=[max(Jaccard_1) min(find(Jaccard_1==max(Jaccard_1))) C_high_threshold_1(min(find(Jaccard_1==max(Jaccard_1))))]';
    
%     dlmwrite(['Euclidean_distance_'  filesname_1_1 '.txt'],euclidean_distance);   
    
    for poly_number=1:9
        poly_number_1=['poly' num2str(poly_number)];   
        for jj=1:size(distance_1,2)            
            % Curve fitting
            % sigmoid (4 parameters)
            cd(codepath);
            [cf_sigmoid_4_parameters,G_sigmoid_4_parameters]=L4P(C_high_threshold_1',distance_1(:,jj));
            cd(filepath);
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         scatter(C_low_threshold_1',distance(:,jj))
    %         hold on
    %         plot(cf_sigmoid_4_parameters)
            features_name={'1_entropy','2_skewness','3_kurtosis','4_autocorrelation','5_cluster_prominence','6_cluster_shade','7_dissimilarity','8_maximum_probability','9_angular_second_moment','10_contrast','11_correlation','12_variance','13_inverse_difference_moment','14_sum_average','15_sum_variance','16_sum_entropy','17_entropy_haralick','18_difference_variance','19_difference_entropy','20_information_measure_of_correlation_1','21_information_measure_of_correlation_2'};
         
            %  saveas(gcf,[filesname_1_1 '_distance_and_fitting_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
            [fx_sigmoid_4_parameters,fxx_sigmoid_4_parameters] = differentiate(cf_sigmoid_4_parameters,C_high_threshold_1);  % 對 object 的微分，cf 是 object，C_low_threshold 是微分範圍，fx 是一次微分的結果，fxx 是二次微分的結果
            fx_sigmoid_4_parameters_1(:,jj)=fx_sigmoid_4_parameters;
            fxx_sigmoid_4_parameters_1(:,jj)=fxx_sigmoid_4_parameters;
            fx_sigmoid_4_parameters=fix(fx_sigmoid_4_parameters);  % 朝 0 的方向取整
            fxx_sigmoid_4_parameters=fix(fxx_sigmoid_4_parameters);  % 朝 0 的方向取整
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fx_sigmoid_4_parameters);
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fxx_sigmoid_4_parameters);
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
            cf_points_sigmoid_4_parameters_0=cf_sigmoid_4_parameters(C_high_threshold_1');
            cf_points_sigmoid_4_parameters(:,jj)=cf_points_sigmoid_4_parameters_0;
            distance_fitting_r2_sigmoid_4_parameters(jj)=G_sigmoid_4_parameters.rsquare;
    %         location_max_fitting_slope_slope_0=max(fxx);
    %         location_max_fitting_slope_slope=round(mean(find(fxx==location_max_fitting_slope_slope_0)));
    %         location_max_fitting_slope_slope_1(jj)=location_max_fitting_slope_slope;
    %         location_max_fitting_slope_slope_2(jj)=C_low_threshold_1(location_max_fitting_slope_slope);
    %         location_max_fitting_slope_slope_3=[location_max_fitting_slope_slope_1;location_max_fitting_slope_slope_2]';
            % For max fitting_slope_slope
    %         location_max_fitting_slope_slope_sigmoid_4_parameters_0=find(fx_sigmoid_4_parameters.*fxx_sigmoid_4_parameters<=0);
    %         max_fitting_slope_slope_points_sigmoid_4_parameters=max(abs(fxx_sigmoid_4_parameters(location_max_fitting_slope_slope_sigmoid_4_parameters_0)));
    %         location_max_fitting_slope_slope_sigmoid_4_parameters_1(jj)=find(abs(fxx_sigmoid_4_parameters)==max_fitting_slope_slope_points_sigmoid_4_parameters);
    %         location_max_fitting_slope_slope_sigmoid_4_parameters_2(jj)=C_low_threshold_1(location_max_fitting_slope_slope_sigmoid_4_parameters_1(jj));
    %         location_max_fitting_slope_slope_sigmoid_4_parameters_3=[location_max_fitting_slope_slope_sigmoid_4_parameters_1;location_max_fitting_slope_slope_sigmoid_4_parameters_2]';
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fxx_sigmoid_4_parameters);
    %         hold on
    %         scatter(location_max_fitting_slope_slope_sigmoid_4_parameters_2(jj),max_fitting_slope_slope_points_sigmoid_4_parameters,120,'*')
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_slope_and_location_max_fitting_slope_slope' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         scatter(C_low_threshold_1',distance(:,jj));
    %         hold on
    %         scatter(location_max_fitting_slope_slope_sigmoid_4_parameters_2(jj),distance(location_max_fitting_slope_slope_sigmoid_4_parameters_1(jj),jj),120,'*')
    %         saveas(gcf,[filesname_1_1 '_distance_and_location_max_fitting_slope_slope' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
            % For min fitting_slope
            location_x_max_fitting_slope_slope_sigmoid_4_parameters_0=find(fx_sigmoid_4_parameters.*fxx_sigmoid_4_parameters<0);
            TF_5_sigmoid_4_parameters_0 = isempty(location_x_max_fitting_slope_slope_sigmoid_4_parameters_0);  %  判斷是否為空集合，是的話，回傳 1，不是，回傳 0 
            if TF_5_sigmoid_4_parameters_0==1
                location_x_max_fitting_slope_slope_sigmoid_4_parameters_0=1;
            else
            end
            max_fitting_slope_slope_slope_points_sigmoid_4_parameters=max(abs(fxx_sigmoid_4_parameters(location_x_max_fitting_slope_slope_sigmoid_4_parameters_0)));
            location_x_max_fitting_slope_slope_sigmoid_4_parameters_1=min(find(abs(fxx_sigmoid_4_parameters)==max_fitting_slope_slope_slope_points_sigmoid_4_parameters));      
            min_fitting_slope_points_sigmoid_4_parameters=min(abs(fx_sigmoid_4_parameters));        
            location_x_min_fitting_slope_sigmoid_4_parameters_0=find(abs(fx_sigmoid_4_parameters)==min_fitting_slope_points_sigmoid_4_parameters);     
            location_x_min_fitting_slope_sigmoid_4_parameters_0=location_x_min_fitting_slope_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_0>location_x_max_fitting_slope_slope_sigmoid_4_parameters_1);
            TF_4_sigmoid_4_parameters_0 = isempty(location_x_min_fitting_slope_sigmoid_4_parameters_0);  %  判斷是否為空集合，是的話，回傳 1，不是，回傳 0        
            if TF_4_sigmoid_4_parameters_0==1
                location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
            else
                location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=min(location_x_min_fitting_slope_sigmoid_4_parameters_0);           
            end
            max_fitting_slope_points_sigmoid_4_parameters=max(abs(fx_sigmoid_4_parameters));
            location_x_max_fitting_slope_sigmoid_4_parameters_0=min(find(abs(fx_sigmoid_4_parameters)==max_fitting_slope_points_sigmoid_4_parameters));        

    %         distance_location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=distance(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj);        
    %         location_min_fitting_slope_sigmoid_4_parameters=[location_x_min_fitting_slope_sigmoid_4_parameters_1;location_x_min_fitting_slope_sigmoid_4_parameters_2;location_y_min_fitting_slope_sigmoid_4_parameters;distance_location_x_min_fitting_slope_sigmoid_4_parameters_1]';
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fx_sigmoid_4_parameters);
    %         hold on
    %         scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),min_fitting_slope_points_sigmoid_4_parameters,120,'*')
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_and_location_x_min_fitting_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fxx_sigmoid_4_parameters);
    %         hold on
    %         scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),min_fitting_slope_points_sigmoid_4_parameters,120,'*')
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_slope_and_location_x_min_fitting_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         scatter(C_low_threshold_1',distance(:,jj));
    %         hold on
    %         plot(C_low_threshold_1',cf_points_sigmoid_4_parameters_0,'r');
    %         if TF_2_sigmoid_4_parameters_0==0
    %              scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),distance(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),120,'*')
    %         else
    %         end
    %         if TF_2_sigmoid_4_parameters_0==0
    %             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),120,'*')
    %         else
    %         end
    %         saveas(gcf,[filesname_1_1 '_distance_and_location_x_min_fitting_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',cf_points_sigmoid_4_parameters_0);
    %         hold on
    %         if TF_2_sigmoid_4_parameters_0==0
    %             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)))
    %         else
    %         end
    %         saveas(gcf,[filesname_1_1 '_fitting_and_location_x_min_fitting_slope_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔        
    %         close all        

            TF_1_sigmoid_4_parameters=isnan(location_x_min_fitting_slope_sigmoid_4_parameters_1);
            if TF_1_sigmoid_4_parameters==1
                 location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
            else
                if location_x_max_fitting_slope_sigmoid_4_parameters_0<=location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)           
                    location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)=location_x_max_fitting_slope_sigmoid_4_parameters_0-abs(location_x_max_fitting_slope_sigmoid_4_parameters_0-location_x_min_fitting_slope_sigmoid_4_parameters_1(jj));
                else
                    location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
                end
            end       
            if location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)<=0
                location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
            else
            end
            TF_3_sigmoid_4_parameters_0=isnan(location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj));        
            if TF_3_sigmoid_4_parameters_0==1
                location_x_mirror_min_fitting_slope_sigmoid_4_parameters_2(jj)=NaN;
            else
                location_x_mirror_min_fitting_slope_sigmoid_4_parameters_2(jj)=C_high_threshold_1(location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj));
            end
            if TF_3_sigmoid_4_parameters_0==1
                location_y_mirror_min_fitting_slope_sigmoid_4_parameters(jj)=NaN;
            else
                location_y_mirror_min_fitting_slope_sigmoid_4_parameters(jj)=cf_points_sigmoid_4_parameters_0(location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)); 
            end
            location_mirror_min_fitting_slope_sigmoid_4_parameters=[location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1;location_x_mirror_min_fitting_slope_sigmoid_4_parameters_2;location_y_mirror_min_fitting_slope_sigmoid_4_parameters]';
    %         location_x_width_step_fitting_slope(jj)=2*abs(location_x_max_fitting_slope_sigmoid_4_parameters_0-location_x_min_fitting_slope_sigmoid_4_parameters_1(jj));
    %         location_x_width_step_fitting_slope=location_x_width_step_fitting_slope';

    %         TF_1_sigmoid_4_parameters_0 = isempty(location_x_min_fitting_slope_sigmoid_4_parameters_0);  %  判斷是否為空集合，是的話，回傳 1，不是，回傳 0
    %         if TF_1_sigmoid_4_parameters_0==1
    %             location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
    %         else           
    %         end
            if TF_3_sigmoid_4_parameters_0==1
                location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
            else
            end        
            TF_2_sigmoid_4_parameters_0=isnan(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj));  %  判斷是否為 NaN，是的話，回傳 1，不是，回傳 0
            if TF_2_sigmoid_4_parameters_0==1
                location_x_min_fitting_slope_sigmoid_4_parameters_2(jj)=NaN;
            else
                location_x_min_fitting_slope_sigmoid_4_parameters_2(jj)=C_high_threshold_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj));
            end
            if TF_2_sigmoid_4_parameters_0==1
                location_y_min_fitting_slope_sigmoid_4_parameters(jj)=NaN;
                distance_location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=NaN;
            else
                location_y_min_fitting_slope_sigmoid_4_parameters(jj)=cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj));
                distance_location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)=distance_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj);
            end
            location_min_fitting_slope_sigmoid_4_parameters=[location_x_min_fitting_slope_sigmoid_4_parameters_1;location_x_min_fitting_slope_sigmoid_4_parameters_2;location_y_min_fitting_slope_sigmoid_4_parameters;distance_location_x_min_fitting_slope_sigmoid_4_parameters_1]';

            figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
            plot(C_high_threshold_1',fx_sigmoid_4_parameters);
            hold on
            if TF_2_sigmoid_4_parameters_0==0
                scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),fx_sigmoid_4_parameters(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
            else
            end
            if TF_3_sigmoid_4_parameters_0==0
                scatter(location_x_mirror_min_fitting_slope_sigmoid_4_parameters_2(jj),fx_sigmoid_4_parameters(location_x_mirror_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
            else
            end
            saveas(gcf,[filesname_1_1 '_fitting_and_locations_min_slope_mirror_' features_name{jj} '_sigmoid_4_parameters.jpg']);  % 圖片存檔        

    %         % Step change detection
    %         for sss=1:size(cf_points_sigmoid_4_parameters_0,1)-1
    %             std_1_sigmoid_4_parameters(sss)=std((cf_points_sigmoid_4_parameters_0(1:sss)));
    %             std_2_sigmoid_4_parameters(sss)=std((cf_points_sigmoid_4_parameters_0(sss:end)));
    %         end
    %         std_3_sigmoid_4_parameters(:,jj)=std_1_sigmoid_4_parameters'+std_2_sigmoid_4_parameters';  % step change 在值最小的位置
    %         
    %         std_3_sigmoid_4_parameters_max=max(std_3_sigmoid_4_parameters);
    %         std_3_sigmoid_4_parameters_min=min(std_3_sigmoid_4_parameters);
    %         for sss=1:size(std_3_sigmoid_4_parameters,2)
    %             figure
    %             plot(std_3_sigmoid_4_parameters(:,sss))
    %             hold on
    %             scatter(find(std_3_sigmoid_4_parameters(:,sss)==std_3_sigmoid_4_parameters_max(sss)),std_3_sigmoid_4_parameters_max(sss))
    %             scatter(find(std_3_sigmoid_4_parameters(:,sss)==std_3_sigmoid_4_parameters_min(sss)),std_3_sigmoid_4_parameters_min(sss))
    %             std_3_sigmoid_4_parameters_slope_max_min(sss)=(std_3_sigmoid_4_parameters_max(sss)-std_3_sigmoid_4_parameters_min(sss))/(find(std_3_sigmoid_4_parameters(:,sss)==std_3_sigmoid_4_parameters_max(sss))-find(std_3_sigmoid_4_parameters(:,sss)==std_3_sigmoid_4_parameters_min(sss)));
    %         end

            % poly
            [xData, yData] = prepareCurveData( C_high_threshold_1', distance_1(:,jj) );
            ft = fittype( poly_number_1 );      
            [cf_poly, G_poly] = fit( xData, yData, ft );        
    %         [cf_exp2,G_exp2] = fit(C_low_threshold_1',distance(:,jj),'exp2');
    %         cd(filepath);
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         scatter(C_low_threshold_1',distance(:,jj))
    %         hold on
    %         plot(cf_exp2)
    %         features_name={'1_intensity','2_entropy','3_skewness','4_kurtosis','5_contrast','6_homogeneity','7_correlation','8_energy','9_angular_second_moment','10_contrast_haralick','11_correlation_haralick','12_variance','13_inverse_difference_moment','14_sum_average','15_sum_variance','16_sum_entropy','17_entropy_haralick','18_difference_variance','19_difference_entropy','20_information_measure_of_correlation_1','21_information_measure_of_correlation_2'};
    %         saveas(gcf,[filesname_1_1 '_distance_and_fitting_' features_name{jj} '_exp2.jpg']);  % 圖片存檔
            [fx_poly,fxx_poly] = differentiate(cf_poly,C_high_threshold_1);  % 對 object 的微分，cf 是 object，C_low_threshold 是微分範圍，fx 是一次微分的結果，fxx 是二次微分的結果
            fx_poly_1(:,jj)=fx_poly;
            fxx_poly_1(:,jj)=fxx_poly;
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fx_exp2);
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_' features_name{jj} '_exp2.jpg']);  % 圖片存檔
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fxx_exp2);
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_slope_' features_name{jj} '_exp2.jpg']);  % 圖片存檔
            cf_points_poly_0=cf_poly(C_high_threshold_1');
            cf_points_poly(:,jj)= cf_points_poly_0;
            distance_fitting_r2_poly(jj)=G_poly.rsquare;
    %         location_max_fitting_slope_slope_exp2_0=find(fx_exp2.*fxx_exp2<=0);
    %         max_fitting_slope_slope_points_exp2=max(abs(fxx_exp2(location_max_fitting_slope_slope_exp2_0)));
    %         location_max_fitting_slope_slope_exp2_1(jj)=find(abs(fxx_exp2)==max_fitting_slope_slope_points_exp2);
    %         location_max_fitting_slope_slope_exp2_2(jj)=C_low_threshold_1(location_max_fitting_slope_slope_exp2_1(jj));
    %         location_max_fitting_slope_slope_exp2_3=[location_max_fitting_slope_slope_exp2_1;location_max_fitting_slope_slope_exp2_2]';
    %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    %         plot(C_low_threshold_1',fxx_exp2);
    %         hold on
    %         scatter(location_max_fitting_slope_slope_exp2_2(jj),max_fitting_slope_slope_points_exp2)
    %         saveas(gcf,[filesname_1_1 '_fitting_slope_slope_and_location_max_fitting_slope_slope' features_name{jj} '_exp2.jpg']);  % 圖片存檔
    %         close all
            % 找原始數據與趨勢線的交點
    %         DD_thrshold=0.001;
    %         DD=find(abs(distance(:,jj)-cf_points_exp2_0)<=DD_thrshold);
    % %         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    % %         scatter(C_low_threshold_1',distance(:,jj))
    % %         hold on
    % %         plot(C_low_threshold_1',cf_points_exp2_0)
    % %         scatter(C_low_threshold_1(DD)',cf_points_exp2_0(DD))
    %         location_min_intersection_cf_fitting_exp2_1(jj)=min(DD);
    %         location_min_intersection_cf_fitting_exp2_2(jj)=C_low_threshold_1(location_min_intersection_cf_fitting_exp2_1(jj));
    %         location_min_intersection_cf_fitting_exp2_3=[location_min_intersection_cf_fitting_exp2_1;location_min_intersection_cf_fitting_exp2_2]';
            L1X = 1:size(distance_1,1); 
            L2X = L1X;
            L1Y=cf_points_poly_0;
            L2Y=distance_1(:,jj);
            cd(codepath)
            [X0,Y0] = intersections(L1X,L1Y,L2X,L2Y);
            cd(filepath)
            X0_greater_than_fx_slope_min=min(round(X0(find(X0>location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)))));
            TF_1_poly_0 = isempty(X0_greater_than_fx_slope_min);  %  判斷是否為空集合，是的話，回傳 1，不是，回傳 0
            if TF_1_poly_0==1
                location_x_appropriate_intersection_distance_fitting_1(jj)=NaN;
            else
                location_x_appropriate_intersection_distance_fitting_1(jj)=X0_greater_than_fx_slope_min;
            end
            if TF_1_poly_0==1
                location_x_appropriate_intersection_distance_fitting_2(jj)=NaN;
            else
                location_x_appropriate_intersection_distance_fitting_2(jj)=C_high_threshold_1(location_x_appropriate_intersection_distance_fitting_1(jj));
            end
            if TF_1_poly_0==1
                location_y_appropriate_intersection_distance_fitting(jj)=NaN;
                distance_location_y_appropriate_intersection_distance_fitting(jj)=NaN;
            else
                location_y_appropriate_intersection_distance_fitting(jj)=cf_points_poly_0(location_x_appropriate_intersection_distance_fitting_1(jj));
                distance_location_y_appropriate_intersection_distance_fitting(jj)=distance_1(location_x_appropriate_intersection_distance_fitting_1(jj),jj);
            end
            location_appropriate_intersection_distance_fitting=[location_x_appropriate_intersection_distance_fitting_1;location_x_appropriate_intersection_distance_fitting_2;location_y_appropriate_intersection_distance_fitting]';

            location_x_optimal_1(jj)=(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)+location_x_appropriate_intersection_distance_fitting_1(jj))/2;
            
            %  去除 outliers
            if outlier_selection==1              
                if location_x_optimal_1(jj)>size(C_high_threshold_1,2)*0.75;
                    location_x_optimal_1(jj)=NaN;
                else
                end
            end
            
            location_x_optimal_round(jj)=round(location_x_optimal_1(jj));
            if location_x_optimal_round(jj)+1>=size(C_high_threshold_1,2)
                location_x_optimal_round(jj)=location_x_optimal_round(jj)-1;
            else
            end
            if location_x_optimal_round(jj)<1
                location_x_optimal_round(jj)=2
            else
            end
            TF_2_poly_0=isnan(location_x_optimal_1(jj));            
            if TF_2_poly_0==1
                location_x_optimal_2(jj)=NaN;
                location_y_optimal(jj)=NaN;
                Precision_2(jj)=NaN;
                meanError_2(jj)=NaN;
                totalCost_2(jj)=NaN;
            else
                if location_x_optimal_1(jj)==location_x_optimal_round(jj)
                    location_x_optimal_2(jj)=C_high_threshold_1(location_x_optimal_1(jj));
                    location_y_optimal(jj)=distance_1(location_x_optimal_1(jj),jj);    
                    Precision_2(jj)=Precision_1(location_x_optimal_1(jj));
                    meanError_2(jj)=meanError_1(location_x_optimal_1(jj));
                    totalCost_2(jj)=totalCost_1(location_x_optimal_1(jj));
                else
                    location_x_optimal_2(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,C_high_threshold_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1),location_x_optimal_1(jj));
                    location_y_optimal(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,distance_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,jj),location_x_optimal_1(jj));
                    Precision_2(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,Precision_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1),location_x_optimal_1(jj));
                    meanError_2(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,meanError_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1),location_x_optimal_1(jj));
                    totalCost_2(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,totalCost_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1),location_x_optimal_1(jj));
                end
            end
            location_optimal=[location_x_optimal_1;location_x_optimal_2;location_y_optimal]';
            if isnan(location_x_optimal_2(jj))==1;
                Dice=NaN;
                Jaccard=NaN;
            else
%                 BW = edge(I,'Canny',[location_x_optimal_2(jj)*ThresholdRatio location_x_optimal_2(jj)], C_sigma);
                cd(codepath);
                BW = canny_detector(I, C_sigma, [location_x_optimal_2(jj)*ThresholdRatio location_x_optimal_2(jj)]);
                cd(filepath);                
                sel1=strel('square',3);
                sel2=strel('square',5);
                BW=imdilate(BW,sel2);
                BW = imfill(BW, 'holes');
                imLabel = bwlabel(BW);% ??通?域?行?
                stats = regionprops(imLabel,'Area');
                [b,index]=sort([stats.Area],'descend');
                if length(stats)<size_order_threshold  % 前幾大
                    bw=imLabel;
                else
                    bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
                end
                cd(codepath);
                [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
                cd(filepath);
            end
            Dice_optimal_individual(poly_number,jj)=Dice;  
            Jaccard_optimal_individual(poly_number,jj)=Jaccard;
            Precision_optimal_individual(poly_number,jj)=Precision_2(jj);
            meanError_optimal_individual(poly_number,jj)=meanError_2(jj);
            totalCost_optimal_individual(poly_number,jj)=totalCost_2(jj);     
           
            figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
            scatter(C_high_threshold_1',distance_1(:,jj))
            hold on
            plot(C_high_threshold_1',cf_points_sigmoid_4_parameters_0,'r')
            if TF_2_sigmoid_4_parameters_0==0
                scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
                scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),distance_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),60,'filled')            
            else
            end        
            plot(C_high_threshold_1',cf_points_poly_0,'k')        
            if TF_1_poly_0==1
            else
                scatter(location_x_appropriate_intersection_distance_fitting_2(jj),cf_points_poly_0(location_x_appropriate_intersection_distance_fitting_1(jj)),60,'filled')        
            end
            TF_1_0=isnan(location_x_optimal_1(jj));
            if TF_1_0==1
            else
            scatter(location_x_optimal_2(jj),location_y_optimal(jj),60,'filled')
            end
            saveas(gcf,[filesname_1_1 '_distance_fitting_locations_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔  
            scatter(Dice_index_max(i,3),distance_1(Dice_index_max(i,2),jj),60,'c','filled')
            saveas(gcf,[filesname_1_1 '_distance_fitting_locations_max_Dice_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔
            distance_1_max_Dice(i,jj)=distance_1(Dice_index_max(i,2),jj);
            close all
            figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
            scatter(C_high_threshold_1',distance_1(:,jj))
            hold on
            plot(C_high_threshold_1',cf_points_sigmoid_4_parameters_0,'r')
            if TF_2_sigmoid_4_parameters_0==0
                scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
                scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),distance_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),60,'filled')            
            else
            end        
            plot(C_high_threshold_1',cf_points_poly_0,'k')        
            if TF_1_poly_0==1
            else
                scatter(location_x_appropriate_intersection_distance_fitting_2(jj),cf_points_poly_0(location_x_appropriate_intersection_distance_fitting_1(jj)),60,'filled')        
            end
            TF_1_0=isnan(location_x_optimal_1(jj));
            if TF_1_0==1
            else
            scatter(location_x_optimal_2(jj),location_y_optimal(jj),60,'filled')
            end
            scatter(Jaccard_index_max(i,3),distance_1(Jaccard_index_max(i,2),jj),60,'c','fileed')
            saveas(gcf,[filesname_1_1 '_distance_fitting_locations_max_Jaccard_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔 
            distance_1_max_Jaccard(i,jj)=distance_1(Jaccard_index_max(i,2),jj);
            close all           
        end        
%         estimation_time_1(poly_number)=toc;
        
%         dlmwrite(['fitting_points_'  filesname_1_1 '_' num2str(C_sigma^2) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '_sigmoid_4_parameters.txt'],cf_points_sigmoid_4_parameters);
        dlmwrite(['fitting_points_'  filesname_1_1 '_sigmoid_4_parameters.txt'],cf_points_sigmoid_4_parameters);
        dlmwrite(['fitting_slope_'  filesname_1_1 '_sigmoid_4_parameters.txt'],fx_sigmoid_4_parameters_1);
        dlmwrite(['fitting_slope_slope_' filesname_1_1 '_sigmoid_4_parameters.txt'],fxx_sigmoid_4_parameters_1); 
    %     dlmwrite(['location_max_curve_slope_slope_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_max_fitting_slope_slope_sigmoid_4_parameters_3); 
    %     dlmwrite(['location_x_min_curve_slope_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_x_min_fitting_slope_sigmoid_4_parameters_3);
        dlmwrite(['fitting_r2_'  filesname_1_1 '_sigmoid_4_parameters.txt'],distance_fitting_r2_sigmoid_4_parameters);
        dlmwrite(['locations_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_min_fitting_slope_sigmoid_4_parameters);
        dlmwrite(['locations_mirror_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_mirror_min_fitting_slope_sigmoid_4_parameters);    

        dlmwrite(['fitting_points_' filesname_1_1 '_' poly_number_1 '.txt'],cf_points_poly);
        dlmwrite(['fitting_slope_'  filesname_1_1 '_'  poly_number_1 '.txt'],fx_poly_1);
        dlmwrite(['fitting_slope_slope_'  filesname_1_1 '_'  poly_number_1 '.txt'],fxx_poly_1); 
    %     dlmwrite(['location_max_curve_slope_slope_'  filesname_1_1 '_' num2str(C_sigma^2) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '_poly4.txt'],location_max_fitting_slope_slope_exp2_3); 
        dlmwrite(['fitting_r2_'  filesname_1_1 '_' poly_number_1 '.txt'],distance_fitting_r2_poly); 
        dlmwrite(['locations_'  filesname_1_1 '_'  poly_number_1 '.txt'],location_appropriate_intersection_distance_fitting); 

        dlmwrite(['location_optimal_'  filesname_1_1 '_' poly_number_1 '.txt'],location_optimal); 
%         dlmwrite(['Estimation_time_fitting_' filesname_1_1 '.txt'],estimation_time_1);
        
        % Performance evaluation for optimal thresholds
        clear location_optimal_0 location_optimal_1 location_optimal_2 thisBoundary thisBoundary_x thisBoundary_y p q Precision Recall F1 TP FP FN assignment meanError totalCost
        location_optimal_0=location_optimal(:,2);
        location_optimal_1=location_optimal_0(~isnan(location_optimal_0));
        location_optimal_2=mean(location_optimal_1);
        if isnan(location_optimal_2)==0
%             BW = edge(I,'Canny',[location_optimal_2*ThresholdRatio location_optimal_2], C_sigma);
            cd(codepath);
            BW = canny_detector(I, C_sigma, [location_optimal_2*ThresholdRatio location_optimal_2]);
            cd(filepath);   
            sel1=strel('square',3);
            sel2=strel('square',5);
            BW=imdilate(BW,sel2);
            BW = imfill(BW, 'holes');
            imLabel = bwlabel(BW);% ??通?域?行?
            stats = regionprops(imLabel,'Area');
            [b,index]=sort([stats.Area],'descend');
            if length(stats)<size_order_threshold  % 前幾大
                bw=imLabel;
            else
                bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
            end
            cd(codepath);
            [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
            cd(filepath);
            Dice_optimal(poly_number)=Dice;  
            Jaccard_optimal(poly_number)=Jaccard;
            bw_edge=bwboundaries(bw);
            k=size(bw_edge,1);
            thisBoundary = bw_edge{k};
            thisBoundary_y=thisBoundary(:,2);
            thisBoundary_x=thisBoundary(:,1);
            figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
            I_1=I_0;
            if I_dimension==3
                if image_size(3)==4;
                    I_1=I_1(:,:,1:3);
                    if dynamic_class_0==2^16;
                        I_1(I_1>=2500)=2500;  % for brightfield image
                        imagesc(I_1)
                         colormap(gray)
                    else
                        imshow(I_1)
                    end
                elseif image_size(3)==3;
                    if dynamic_class_0==2^16;
                        I_1(I_1>=2500)=2500;  % for brightfield image
                        imagesc(I_1)
                         colormap(gray)
                    else
                        imshow(I_1)
                    end
                end
            else
                if dynamic_class_0==2^16;
                    I_1(I_1>=2500)=2500;  % for brightfield image
                    imagesc(I_1)
                    colormap(gray)
                else
                    imshow(I_1)
                end
            end
            axis image
            axis off
            hold on
            plot(thisBoundary_y, thisBoundary_x,'r:', 'LineWidth', 2)
            saveas(gcf,[filesname_1_1 '_optimal_' poly_number_1 '_with_ROI.jpg']);  % 圖片存檔
            
            % Performance evaluation of edge edtection
            p=load(['ROI_' filesname_1_1 '.txt']);
            q=[thisBoundary_y thisBoundary_x];
            cd(codepath);
            p = resampleBoundary(p, 200);
            q = resampleBoundary(q, 200);        
            [Precision, Recall, F1, TP, FP, FN, assignment, meanError, totalCost] = evaluateBoundaryMunkres(p, q, size(I), d_ratio, mode);
            cd(filepath)
            Precision_optimal(poly_number)=Precision;
            meanError_optimal(poly_number)=meanError;
            totalCost_optimal(poly_number)=totalCost;
        else
        end
        
    %     % PCA
    %     parameters=[ROI_parameters Background_parameters];
    %     parameters=zscore(parameters);  % 標準化(Z 分數)
    %     parameters=parameters';
    %     [coeff,score,latent,tsquared,explained,mu]=princomp(parameters);
    %     %         parameters_0 = bsxfun(@minus,parameters_0,mean(parameters_0));  % 對兩個矩陣A和B之間的每一個元素進行指定的計算
    %     %         parameters_0_std=std(parameters_0);
    %     %         for ii=1:size(parameters_0,2)
    %     %             parameters(:,ii)=parameters_0(:,ii)/parameters_0_std(1,ii);
    %     %         end
    %     dlmwrite(['PCA_score_'  filesname_1_1 '.txt'],score);
    %     dlmwrite(['PCA_explained_'  filesname_1_1 '.txt'],explained);
    end
    Jaccard_index_optimal(:,i)=Jaccard_optimal';
    Dice_index_optimal(:,i)=Dice_optimal';
    Precision_index_optimal(:,i)=Precision_optimal';
    meanError_index_optimal(:,i)=meanError_optimal';
    totalCost_index_optimal(:,i)=totalCost_optimal';
    dlmwrite(['Dice_index_optimal_individual_feature_'  filesname_1_1 '.txt'],Dice_optimal_individual);
    dlmwrite(['Jaccard_index_optimal_individual_feature_'  filesname_1_1 '.txt'],Jaccard_optimal_individual);
    dlmwrite(['Precision_optimal_individual_feature_'  filesname_1_1 '.txt'],Precision_optimal_individual);
    dlmwrite(['totalCost_optimal_individual_feature_'  filesname_1_1 '.txt'],totalCost_optimal_individual);
    dlmwrite(['meanError_optimal_individual_feature_'  filesname_1_1 '.txt'],meanError_optimal_individual);
     
    % Performance evaluation for optimal average thresholds
    clear poly_number_1 location_optimal_0 location_optimal_1 location_optimal_2 location_optimal_average_location_3_x location_optimal_average_location_3_y thisBoundary thisBoundary_x thisBoundary_y 
    for nn=1:9
        poly_number_1=['poly' num2str(nn)];
        location_optimal_average_0=load(['location_optimal_'  filesname_1_1 '_' poly_number_1 '.txt']);
        location_optimal_average_1(nn,:)=location_optimal_average_0(:,2)';
%         location_optimal_average_location_1(nn,:)=location_optimal_average_0(:,1)';
    end
    for nnn=1:size(location_optimal_average_1,2)
        location_optimal_average_2=location_optimal_average_1(:,nnn);
        location_optimal_average_2=location_optimal_average_2(~isnan(location_optimal_average_2));
        location_optimal_average_3(nnn)=mean(location_optimal_average_2);
%         location_optimal_average_location_2=location_optimal_average_location_1(:,nnn);
%         location_optimal_average_location_2=location_optimal_average_location_2(~isnan(location_optimal_average_location_2));
%         location_optimal_average_location_3(nnn)=mean(location_optimal_average_location_2);
    end    
    for nnnn=1:size(location_optimal_average_1,2)
%         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%         scatter(C_high_threshold_1',distance_1(:,nnnn))
%         hold on
%         if isnan(location_optimal_average_location_3(nnnn))==1
%             location_optimal_average_location_3_x(nnnn)=NaN;
%             location_optimal_average_location_3_y(nnnn)=NaN;
%         else
%             location_optimal_average_location_3_x(nnnn)=interp1(round(location_optimal_average_location_3(nnnn))-1:round(location_optimal_average_location_3(nnnn))+1,C_high_threshold_1(round(location_optimal_average_location_3(nnnn))-1:round(location_optimal_average_location_3(nnnn))+1),location_optimal_average_location_3(nnnn));
%             location_optimal_average_location_3_y(nnnn)=interp1(round(location_optimal_average_location_3(nnnn))-1:round(location_optimal_average_location_3(nnnn))+1,distance_1(round(location_optimal_average_location_3(nnnn))-1:round(location_optimal_average_location_3(nnnn))+1),location_optimal_average_location_3(nnnn));
%             scatter(location_optimal_average_location_3_x(nnnn),location_optimal_average_location_3_y(nnnn),60,'r','filled')
%         end
%         saveas(gcf,[filesname_1_1 '_distance_fitting_location_optimal_average_location_' features_name{nnnn} '.jpg']);  % 圖片存檔
        TF_2_optimal_average=isnan(location_optimal_average_3(nnnn));
        if TF_2_optimal_average==1
            Dice_optimal_average(nnnn,i)=NaN;
            Jaccard_optimal_average(nnnn,i)=NaN;            
        else
%             BW = edge(I,'Canny',[location_optimal_average_3(nnnn)*ThresholdRatio location_optimal_average_3(nnnn)], C_sigma);
            cd(codepath);
            BW = canny_detector(I, C_sigma, [location_optimal_average_3(nnnn)*ThresholdRatio location_optimal_average_3(nnnn)]);
            cd(filepath); 
            sel1=strel('square',3);
            sel2=strel('square',5);
            BW=imdilate(BW,sel2);
            BW = imfill(BW, 'holes');
            imLabel = bwlabel(BW);% ??通?域?行?
            stats = regionprops(imLabel,'Area');
            [b,index]=sort([stats.Area],'descend');
            if length(stats)<size_order_threshold  % 前幾大
                bw=imLabel;
            else
                bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
            end
            cd(codepath);
            [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
            cd(filepath);
            Dice_optimal_average(nnnn,i)=Dice;  
            Jaccard_optimal_average(nnnn,i)=Jaccard;
            bw_edge=bwboundaries(bw);
            k=size(bw_edge,1);
            thisBoundary = bw_edge{k};
            thisBoundary_y=thisBoundary(:,2);
            thisBoundary_x=thisBoundary(:,1);
            figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
            I_1=I_0;
            if I_dimension==3
                if image_size(3)==4;
                    I_1=I_1(:,:,1:3);
                    if dynamic_class_0==2^16;
                        I_1(I_1>=2500)=2500;  % for brightfield image
                        imagesc(I_1)
                        colormap(gray)
                    else
                        imshow(I_1)
                    end
                elseif image_size(3)==3;
                    if dynamic_class_0==2^16;
                        I_1(I_1>=2500)=2500;  % for brightfield image
                        imagesc(I_1)
                        colormap(gray)
                    else
                        imshow(I_1)
                    end
                end
            else
                if dynamic_class_0==2^16;
                    I_1(I_1>=2500)=2500;  % for brightfield image
                    imagesc(I_1)
                    colormap(gray)
                else
                    imshow(I_1)
                end
            end
            axis image
            axis off
            hold on
            plot(thisBoundary_y, thisBoundary_x,'r:', 'LineWidth', 2)
            saveas(gcf,[filesname_1_1 '_optimal_poly_average_' features_name{nnnn} '_with_ROI.jpg']);  % 圖片存檔
        end
    end
%     dlmwrite(['location_optimal_average_location_'  filesname_1_1 '.txt'],[location_optimal_average_location_3_x;location_optimal_average_location_3_y]');
    close all
    
%     location_optimal_average_location_4=location_optimal_average_location_1;
%     location_optimal_average_location_4=location_optimal_average_location_4(~isnan(location_optimal_average_location_4));
%     location_optimal_average_location_5(i)=mean(location_optimal_average_location_4);
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%     distance_2=mean(distance_1');
%     scatter(C_high_threshold_1',mean(distance_1'))
%     hold on
%     if isnan(location_optimal_average_location_5(i))==1
%         location_optimal_average_location_5_x(i)=NaN;
%         location_optimal_average_location_5_y(i)=NaN;
%     else
%         location_optimal_average_location_5_x(i)=interp1(round(location_optimal_average_location_5(i))-1:round(location_optimal_average_location_5(i))+1,C_high_threshold_1(round(location_optimal_average_location_5(i))-1:round(location_optimal_average_location_5(i))+1),location_optimal_average_location_5(i));
%         location_optimal_average_location_5_y(i)=interp1(round(location_optimal_average_location_5(i))-1:round(location_optimal_average_location_5(i))+1,distance_2(round(location_optimal_average_location_5(i))-1:round(location_optimal_average_location_5(i))+1),location_optimal_average_location_5(i));
%         scatter(location_optimal_average_location_5_x(i),location_optimal_average_location_5_y(i),60,'r','filled')
%     end
%     saveas(gcf,[filesname_1_1 'distance_fitting_location_optimal_all_average_location.jpg']);  % 圖片存檔   
    
    location_optimal_average_4=location_optimal_average_1;
    location_optimal_average_4=location_optimal_average_4(~isnan(location_optimal_average_4));
    location_optimal_average_5(i)=mean(location_optimal_average_4);
%     BW = edge(I,'Canny',[location_optimal_average_5(i)*ThresholdRatio location_optimal_average_5(i)], C_sigma);
    cd(codepath);
    BW = canny_detector(I, C_sigma, [location_optimal_average_5(i)*ThresholdRatio location_optimal_average_5(i)]);
    cd(filepath);
    sel1=strel('square',3);
    sel2=strel('square',5);
    BW=imdilate(BW,sel2);
    BW = imfill(BW, 'holes');
    imLabel = bwlabel(BW);% ??通?域?行?
    stats = regionprops(imLabel,'Area');
    [b,index]=sort([stats.Area],'descend');
    if length(stats)<size_order_threshold  % 前幾大
        bw=imLabel;
    else
        bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
    end
    cd(codepath);
    [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
    cd(filepath);
    Dice_optimal_all_average(i)=Dice;  
    Jaccard_optimal_all_average(i)=Jaccard;
    bw_edge=bwboundaries(bw);
    k=size(bw_edge,1);
    thisBoundary = bw_edge{k};
    thisBoundary_y=thisBoundary(:,2);
    thisBoundary_x=thisBoundary(:,1);
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    I_1=I_0;
    if I_dimension==3
        if image_size(3)==4;
            I_1=I_1(:,:,1:3);
            if dynamic_class_0==2^16;
                I_1(I_1>=2500)=2500;  % for brightfield image
                imagesc(I_1)
                colormap(gray)
            else
                imshow(I_1)
            end
        elseif image_size(3)==3;
            if dynamic_class_0==2^16;
                I_1(I_1>=2500)=2500;  % for brightfield image
                imagesc(I_1)
                colormap(gray)
            else
                imshow(I_1)
            end
        end
    else
        if dynamic_class_0==2^16;
            I_1(I_1>=2500)=2500;  % for brightfield image
            imagesc(I_1)
            colormap(gray)
        else
            imshow(I_1)
        end
    end
    axis image
    axis off
    hold on
    plot(thisBoundary_y, thisBoundary_x,'r:', 'LineWidth', 2)
    saveas(gcf,[filesname_1_1 '_optimal_average_all_with_ROI.jpg']);  % 圖片存檔  
    close all
    
    % Performance evaluation for Matlab default
    clear thisBoundary thisBoundary_x thisBoundary_y
    [BW,thres] = edge(I,'Canny');
%     thres_1(i)=thres(2);
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%     scatter(C_high_threshold_1',distance_1(:,jj))
%     hold on
%     [~, iMin] = min(abs(C_high_threshold_1-location_optimal_Matlab_default_x(i)));  % 找最接近值
%     location_optimal_Matlab_default_x(i)=interp1(iMin-1:iMin+1,C_high_threshold_1(iMin-1:iMin+1),C_high_threshold_1(thres(2)));  % 未完成  
%     location_optimal_Matlab_default_y(i)=interp1(iMin-1:iMin+1,distance_2(iMin-1:iMin+1,jj),location_optimal_average_location_5(i));  % 未完成
%     scatter(location_optimal_Matlab_default_x(i),location_optimal_Matlab_default_y(i),60,'r','filled')
%     saveas(gcf,[filesname_1_1 'distance_fitting_location_Matlab_default.jpg']);  % 圖片存檔    
    sel1=strel('square',3);
    sel2=strel('square',5);
    BW=imdilate(BW,sel2);
    BW = imfill(BW, 'holes');
    imLabel = bwlabel(BW);% ??通?域?行?
    stats = regionprops(imLabel,'Area');
    [b,index]=sort([stats.Area],'descend');
    if length(stats)<size_order_threshold  % 前幾大
        bw=imLabel;
    else
        bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
    end
    cd(codepath);
    [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
    cd(filepath);
    Dice_Matlab_default(i)=Dice;  
    Jaccard_Matlab_default(i)=Jaccard; 
    bw_edge=bwboundaries(bw);
    k=size(bw_edge,1);
    thisBoundary = bw_edge{k};
    thisBoundary_y=thisBoundary(:,2);
    thisBoundary_x=thisBoundary(:,1);
    figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
    I_1=I_0;
    if I_dimension==3
        if image_size(3)==4;
            I_1=I_1(:,:,1:3);
            if dynamic_class_0==2^16;
                I_1(I_1>=2500)=2500;  % for brightfield image
                imagesc(I_1)
                colormap(gray)
            else
                imshow(I_1)
            end
        elseif image_size(3)==3;
            if dynamic_class_0==2^16;
                I_1(I_1>=2500)=2500;  % for brightfield image
                imagesc(I_1)
                colormap(gray)
            else
                imshow(I_1)
            end
        end
    else
        if dynamic_class_0==2^16;
            I_1(I_1>=2500)=2500;  % for brightfield image
            imagesc(I_1)
            colormap(gray)
        else
            imshow(I_1)
        end
    end
    axis image
    axis off
    hold on
    plot(thisBoundary_y, thisBoundary_x,'r:', 'LineWidth', 2)
    saveas(gcf,[filesname_1_1 '_Matlab_default' '_with_ROI.jpg']);  % 圖片存檔 
    close all
    estimation_time(i)=toc; 
end

dlmwrite(['Estimation_time.txt'],estimation_time);
dlmwrite(['Entropy_of_whole_image_for_threshold_range_selection.txt'],hs_sum_counts);
% dlmwrite(['Contrast_of_whole_image_for_threshold_range_selection'],[AAA]);
dlmwrite(['location_Dice_index_max.txt'],[Dice_index_max(:,2:3) distance_1_max_Dice]);
dlmwrite(['location_Jaccard_index_max.txt'],[Jaccard_index_max(:,2:3) distance_1_max_Jaccard]);
dlmwrite(['Dice_Jaccard_index_max.txt'],[Dice_index_max(:,1) Jaccard_index_max(:,1)]);
dlmwrite(['Dice_index_Matlab_default.txt'],Dice_Matlab_default);
dlmwrite(['Jaccard_index_Matlab_default.txt'],Jaccard_Matlab_default); 
dlmwrite(['Dice_index_location_optimal.txt'],Dice_index_optimal);
dlmwrite(['Jaccard_index_location_optimal.txt'],Jaccard_index_optimal); 
dlmwrite(['Dice_index_optimal_poly_average.txt'],Dice_optimal_average);
dlmwrite(['Precision_index_location_optimal.txt'],Precision_index_optimal); 
dlmwrite(['meanError_index_location_optimal.txt'],meanError_index_optimal); 
dlmwrite(['totalCost_index_location_optimal.txt'],totalCost_index_optimal); 
dlmwrite(['Jaccard_index_optimal_poly_average.txt'],Jaccard_optimal_average);
dlmwrite(['Dice_index_optimal_all_average.txt'],Dice_optimal_all_average);
dlmwrite(['Jaccard_index_optimal_all_average.txt'],Jaccard_optimal_all_average);
% dlmwrite(['location_optimal_all_average_location_'  filesname_1_1 '.txt'],[location_optimal_average_location_5_x;location_optimal_average_location_5_y]');
% dlmwrite(['location_Matlab_default_'  filesname_1_1 '.txt'],[location_optimal_Matlab_default_x;location_optimal_Matlab_default_y]');

% for poly_number_2=1:9
%     for iii=1:lengthFiles
%         cell_str_2 = strsplit(filesname_1_0(iii).name,'.');  % 分解字串，以 . 為界線
%         filesname_1_2=cell_str_2{1};
%         Dice_integration_1=load(['Dice_index_optimal_individual_feature_'  filesname_1_2 '.txt']);
%         Dice_integration_2(iii,:)=Dice_integration_1(poly_number_2,:);
%         Jaccard_integration_1=load(['Jaccard_index_optimal_individual_feature_'  filesname_1_2 '.txt']);
%         Jaccard_integration_2(iii,:)=Jaccard_integration_1(poly_number_2,:);  
%     end
%     dlmwrite(['Dice_index_integration_poly' num2str(poly_number_2) '.txt'],Dice_integration_2);
%     dlmwrite(['Jaccard_index_integration_poly' num2str(poly_number_2) '.txt'],Jaccard_integration_2);
% end
