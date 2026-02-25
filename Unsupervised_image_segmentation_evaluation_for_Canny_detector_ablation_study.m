clear all
clc
close all
warning off

filepath='C:\Users\User\Desktop\REVISE 2\Ablation study\Replacing texture features with simpler features\Type_XII';
codepath='C:\Users\User\Desktop\程式';
files = dir(fullfile(filepath,'*.tif'));
if isempty(files)==1
    files = dir(fullfile(filepath,'*.png'));
end
if isempty(files)==1
    files = dir(fullfile(filepath,'*.tiff'));
end
lengthFiles = length(files);
% pause(300);

C_sigma=sqrt(1);
ThresholdRatio=0.476;
poly_number_0=6;

C_threshold_number=2^8;
C_threshold_number_threshold_0=0.95;
outlier_selection=1; % 去除 outluers: 1，不去除 outliers: 2
C_threshold_number_threshold=1-C_threshold_number_threshold_0;
size_order_threshold=1;

for i=1:lengthFiles;
    % 每張圖的 threshold range 不同，需要清除相關參數
    clear A AA threshold_range_0 C_high_threshold C_high_threshold_1 distance_1 ROI_parameters Background_parameters Dice_1 Jaccard_1 Dice_Jaccard_index Accuracy Sensitivity Fmeasure Precision MCC Dice Jaccard Specitivity cf_points_sigmoid_4_parameters fx_sigmoid_4_parameters_1 fxx_sigmoid_4_parameters_1 distance_fitting_r2_sigmoid_4_parameters location_min_fitting_slope_sigmoid_4_parameters location_mirror_min_fitting_slope_sigmoid_4_parameters cf_points_poly fx_poly_1 fxx_poly_1 distance_fitting_r2_poly location_appropriate_intersection_distance_fitting location_optimal
    tic;
    kk=1;
    cd(filepath);
    dlmwrite(['000_ThresholdRatio_' num2str(ThresholdRatio) '.txt'], ThresholdRatio);
    dlmwrite(['000_PolyNumber_' num2str(poly_number_0) '.txt'], poly_number_0);
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
    threshold_high_0=threshold_range_0;
    threshold_high_0(threshold_range_0<max(find(counts_2==max(counts_2))))=[];
    threshold_high_0=min(threshold_high_0);
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
    plot(C_high_threshold_0,counts_1)
    hold on
    ylim([min(counts_1(:))-0.05*max(counts_1(:)) max(counts_1(:))+0.05*max(counts_1(:))])
    xlim([-0.025 1.025])
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
    
    dlmwrite([filesname_1_1 '_threshold_range.txt'],[threshold_high_0 C_high_threshold_high threshold_low_0 C_high_threshold_low]);
    dlmwrite(['C_high_threshold_'  filesname_1_1 '_' num2str(C_sigma^2) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '.txt'],C_high_threshold_1');
    
    % 輸入手動分割的結果
    bw_manual_segmentation=load(['RANGE_' filesname_1_1 '.txt']);
    
    for C_high_threshold=C_high_threshold_1
        C_low_threshold=C_high_threshold*ThresholdRatio;
        filesname_1=['_' num2str(C_sigma^2) '_' num2str(C_high_threshold) '_' num2str(C_low_threshold)];
        
        BW = edge(I,'Canny',[C_low_threshold C_high_threshold], C_sigma);        
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
        saveas(gcf,[filesname_1_1 '_1.jpg']);  % 圖片存檔
        hold on
        plot(thisBoundary_y, thisBoundary_x,'b', 'LineWidth', 3)
        saveas(gcf,[filesname_1_1 '_' num2str(kk) filesname_1 '_with_ROI.jpg']);  % 圖片存檔
        
        % ROI parameter
%         % Texture parameter
        I_ROI=I.*logical(bw);
        ROI_intensity_0=I_ROI;
        ROI_intensity_0(ROI_intensity_0==0)=[];
        % 標準化(Z 分數)
        ROI_intensity_0=zscore(ROI_intensity_0);
%         % Shannon entropy
%         histnumber=32;
%         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%         n_ROI=histogram(ROI_intensity_0,histnumber);
%         n_2_ROI=get(n_ROI,'Values');  % 讀取 histogram 相關數據時不能關圖
%         n_1_ROI=n_2_ROI/size(ROI_intensity_0,2); 
%         hs_ROI=(-1)*n_1_ROI.*(log2(n_1_ROI));
%         hs_ROI=hs_ROI(~isnan(hs_ROI));  % 每個長條 in 直方圖
%         hs_sum_ROI=sum(hs_ROI);   % 每個直方圖的範圍
%         % Skwness
%         s_ROI=skewness(ROI_intensity_0);
%         % kurtosis
%         k_texture_ROI=kurtosis(ROI_intensity_0);
%         % GLCM
%         glcms_ROI = graycomatrix(ROI_intensity_0);
%         cd(codepath);
%         [Haralick_ROI] = haralickTextureFeatures(glcms_ROI);
%         [out_ROI] = GLCM_Features1(glcms_ROI);
%         cd(filepath);
%         ROI_parameters_00=[hs_sum_ROI;s_ROI;k_texture_ROI];
%         ROI_parameters_0=[ROI_parameters_00;out_ROI.autoc;out_ROI.cprom;out_ROI.cshad;out_ROI.dissi;out_ROI.maxpr;Haralick_ROI(1:13)];
%         ROI_parameters(:,kk)=ROI_parameters_0;
        ROI_parameters_0=[mean(ROI_intensity_0);std(ROI_intensity_0);min(ROI_intensity_0);median(ROI_intensity_0);max(ROI_intensity_0)];
        ROI_parameters(:,kk)=ROI_parameters_0;
        
        % Background parameter
%          % Texture parameter
        I_background=I.*logical(imcomplement(bw));
        Background_intensity_0=I_background;
        Background_intensity_0(Background_intensity_0==0)=[];        
        % 標準化(Z 分數)
        Background_intensity_0=zscore(Background_intensity_0);        
%         % Shannon entropy
%         histnumber=32;
%         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%         n_background=histogram(Background_intensity_0,histnumber);
%         n_2_background=get(n_background,'Values');  % 讀取 histogram 相關數據時不能關圖
%         n_1_background=n_2_background/size(Background_intensity_0,2); 
%         hs_background=(-1)*n_1_background.*(log2(n_1_background));
%         hs_background=hs_background(~isnan(hs_background));  % 每個長條 in 直方圖
%         hs_sum_background=sum(hs_background);   % 每個直方圖的範圍
%         % Skwness
%         s_background=skewness(Background_intensity_0);
%         % kurtosis
%         k_texture_background=kurtosis(Background_intensity_0);
%         % GLCM
%         glcms_Background = graycomatrix(Background_intensity_0);
%         cd(codepath);
%         [Haralick_Background] = haralickTextureFeatures(glcms_Background);
%         [out_Background] = GLCM_Features1(glcms_Background);
%         cd(filepath);
%         Background_parameters_00=[hs_sum_background;s_background;k_texture_background];
%         Background_parameters_0=[Background_parameters_00;out_Background.autoc;out_Background.cprom;out_Background.cshad;out_Background.dissi;out_Background.maxpr;Haralick_Background(1:13)];
%         Background_parameters(:,kk)=Background_parameters_0;
        Background_parameters_0=[mean(Background_intensity_0);std(Background_intensity_0);min(Background_intensity_0);median(Background_intensity_0);max(Background_intensity_0)];
        Background_parameters(:,kk)=Background_parameters_0;
        
        % Performance evaluation of edge edtection
        cd(codepath);
        [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores(bw_manual_segmentation, bw);
        cd(filepath);
        Dice_1(kk)=Dice;
        Jaccard_1(kk)=Jaccard;
        
        kk=kk+1;
        close all        
    end
    
    distance_1=sqrt((ROI_parameters'-Background_parameters').^2);
    dlmwrite(['ROI_parameters_'  filesname_1_1 '.txt'],ROI_parameters');
    dlmwrite(['Background_parameters_'  filesname_1_1 '.txt'],Background_parameters');
    dlmwrite(['Distance_'  filesname_1_1 '.txt'],distance_1);
    dlmwrite(['Jaccard_index_'  filesname_1_1 '.txt'],Jaccard_1);
    dlmwrite(['Dice_index_'  filesname_1_1 '.txt'],Dice_1);
    
    Dice_index_max(i,:)=[max(Dice_1) min(find(Dice_1==max(Dice_1))) C_high_threshold_1(min(find(Dice_1==max(Dice_1))))]';
    Jaccard_index_max(i,:)=[max(Jaccard_1) min(find(Jaccard_1==max(Jaccard_1))) C_high_threshold_1(min(find(Jaccard_1==max(Jaccard_1))))]';
    
    % Fittings and optimal threshold
    poly_number=poly_number_0;
    poly_number_1=['poly' num2str(poly_number)];
    for jj=1:size(distance_1,2) 
        % Curve fitting
        % sigmoid (4 parameters)
        cd(codepath);
        [cf_sigmoid_4_parameters,G_sigmoid_4_parameters]=L4P(C_high_threshold_1',distance_1(:,jj));
        cd(filepath);
%         features_name={'1_entropy','2_skewness','3_kurtosis','4_autocorrelation','5_cluster_prominence','6_cluster_shade','7_dissimilarity','8_maximum_probability','9_angular_second_moment','10_contrast','11_correlation','12_variance','13_inverse_difference_moment','14_sum_average','15_sum_variance','16_sum_entropy','17_entropy_haralick','18_difference_variance','19_difference_entropy','20_information_measure_of_correlation_1','21_information_measure_of_correlation_2'};
        features_name={'1_mean','2_Std','3_Min','4_median','5_max'};

        
        [fx_sigmoid_4_parameters,fxx_sigmoid_4_parameters] = differentiate(cf_sigmoid_4_parameters,C_high_threshold_1);  % 對 object 的微分，cf 是 object，C_low_threshold 是微分範圍，fx 是一次微分的結果，fxx 是二次微分的結果
        fx_sigmoid_4_parameters_1(:,jj)=fx_sigmoid_4_parameters;
        fxx_sigmoid_4_parameters_1(:,jj)=fxx_sigmoid_4_parameters;
        fx_sigmoid_4_parameters=fix(fx_sigmoid_4_parameters);  % 朝 0 的方向取整
        fxx_sigmoid_4_parameters=fix(fxx_sigmoid_4_parameters);  % 朝 0 的方向取整
        
        cf_points_sigmoid_4_parameters_0=cf_sigmoid_4_parameters(C_high_threshold_1');
        cf_points_sigmoid_4_parameters(:,jj)=cf_points_sigmoid_4_parameters_0;
        distance_fitting_r2_sigmoid_4_parameters(jj)=G_sigmoid_4_parameters.rsquare;
        
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
        
        % Poly fitting
        [xData, yData] = prepareCurveData( C_high_threshold_1', distance_1(:,jj) );
        ft = fittype(poly_number_1);
        [cf_poly, G_poly] = fit(xData, yData, ft); 
        
        [fx_poly,fxx_poly] = differentiate(cf_poly,C_high_threshold_1);  % 對 object 的微分，cf 是 object，C_low_threshold 是微分範圍，fx 是一次微分的結果，fxx 是二次微分的結果
        fx_poly_1(:,jj)=fx_poly;
        fxx_poly_1(:,jj)=fxx_poly;
        
        cf_points_poly_0=cf_poly(C_high_threshold_1');
        cf_points_poly(:,jj)= cf_points_poly_0;
        distance_fitting_r2_poly(jj)=G_poly.rsquare;
        
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
        
        % Optimal threshold
        location_x_optimal_1(jj)=(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)+location_x_appropriate_intersection_distance_fitting_1(jj))/2;
        %  去除 outliers
        if outlier_selection==1
            if location_x_optimal_1(jj)>size(C_high_threshold_1,2)*0.75;
                location_x_optimal_1(jj)=NaN;
            else
            end
        end
        % 異常值處理
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
        else
            if location_x_optimal_1(jj)==location_x_optimal_round(jj)
                location_x_optimal_2(jj)=C_high_threshold_1(location_x_optimal_1(jj));
                location_y_optimal(jj)=distance_1(location_x_optimal_1(jj),jj);
            else
                location_x_optimal_2(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,C_high_threshold_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1),location_x_optimal_1(jj));
                location_y_optimal(jj)=interp1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,distance_1(location_x_optimal_round(jj)-1:location_x_optimal_round(jj)+1,jj),location_x_optimal_1(jj));
            end
        end
        location_optimal=[location_x_optimal_1;location_x_optimal_2;location_y_optimal]';
        
%         if isnan(location_x_optimal_2(jj))==1;
%             Dice=NaN;
%             Jaccard=NaN;
%         else
%             BW = edge(I,'Canny',[location_x_optimal_2(jj)*ThresholdRatio location_x_optimal_2(jj)], C_sigma);
%             sel1=strel('square',3);
%             sel2=strel('square',5);
%             BW=imdilate(BW,sel2);
%             BW = imfill(BW, 'holes');
%             imLabel = bwlabel(BW);% ??通?域?行?
%             stats = regionprops(imLabel,'Area');
%             [b,index]=sort([stats.Area],'descend');
%             if length(stats)<size_order_threshold  % 前幾大
%                 bw=imLabel;
%             else
%                 bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
%             end       
%             cd(codepath);
%             [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores(bw_manual_segmentation,bw);
%             cd(filepath);       
%         end
%         Dice_optimal_individual(jj)=Dice;
%         Jaccard_optimal_individual(jj)=Jaccard;
%         
%         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%         scatter(C_high_threshold_1',distance_1(:,jj))
%         hold on
%         plot(C_high_threshold_1',cf_points_sigmoid_4_parameters_0,'r')
%         if TF_2_sigmoid_4_parameters_0==0
%             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
%             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),distance_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),60,'filled')
%         else
%         end
%         plot(C_high_threshold_1',cf_points_poly_0,'k')
%         if TF_1_poly_0==1
%         else
%             scatter(location_x_appropriate_intersection_distance_fitting_2(jj),cf_points_poly_0(location_x_appropriate_intersection_distance_fitting_1(jj)),60,'filled')
%         end
%         TF_1_0=isnan(location_x_optimal_1(jj));
%         if TF_1_0==1
%         else
%             scatter(location_x_optimal_2(jj),location_y_optimal(jj),60,'filled')
%         end
%         saveas(gcf,[filesname_1_1 '_distance_fitting_locations_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔
%         scatter(Dice_index_max(i,3),distance_1(Dice_index_max(i,2),jj),60,'c','filled')
%         saveas(gcf,[filesname_1_1 '_distance_fitting_locations_max_Dice_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔
%         distance_1_max_Dice(i,jj)=distance_1(Dice_index_max(i,2),jj);
%         close all
%         figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%         scatter(C_high_threshold_1',distance_1(:,jj))
%         hold on
%         plot(C_high_threshold_1',cf_points_sigmoid_4_parameters_0,'r')
%         if TF_2_sigmoid_4_parameters_0==0
%             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),cf_points_sigmoid_4_parameters_0(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj)),60,'filled')
%             scatter(location_x_min_fitting_slope_sigmoid_4_parameters_2(jj),distance_1(location_x_min_fitting_slope_sigmoid_4_parameters_1(jj),jj),60,'filled')
%         else
%         end
%         plot(C_high_threshold_1',cf_points_poly_0,'k')
%         if TF_1_poly_0==1
%         else
%             scatter(location_x_appropriate_intersection_distance_fitting_2(jj),cf_points_poly_0(location_x_appropriate_intersection_distance_fitting_1(jj)),60,'filled')
%         end
%         TF_1_0=isnan(location_x_optimal_1(jj));
%         if TF_1_0==1
%         else
%             scatter(location_x_optimal_2(jj),location_y_optimal(jj),60,'filled')
%         end
%         scatter(Jaccard_index_max(i,3),distance_1(Jaccard_index_max(i,2),jj),60,'c','fileed')
%         saveas(gcf,[filesname_1_1 '_distance_fitting_locations_max_Jaccard_' features_name{jj} '_' poly_number_1 '.jpg']);  % 圖片存檔
%         distance_1_max_Jaccard(i,jj)=distance_1(Jaccard_index_max(i,2),jj);
%         close all         
    end
    
%     location_optimal_all_mean=mean(location_x_optimal_2(~isnan(location_x_optimal_2)));
%     % Imaging using optimal parameters
%     BW = edge(I,'Canny',[location_optimal_all_mean*ThresholdRatio location_optimal_all_mean], C_sigma);
%     sel1=strel('square',3);
%     sel2=strel('square',5);
%     BW=imdilate(BW,sel2);
%     BW = imfill(BW, 'holes');
%     imLabel = bwlabel(BW);% ??通?域?行?
%     stats = regionprops(imLabel,'Area');
%     [b,index]=sort([stats.Area],'descend');
%     if length(stats)<size_order_threshold  % 前幾大
%         bw=imLabel;
%     else
%         bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
%     end
%     bw_edge=bwboundaries(bw);
%     k=size(bw_edge,1);
%     thisBoundary = bw_edge{k};
%     thisBoundary_y=thisBoundary(:,2);
%     thisBoundary_x=thisBoundary(:,1);
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%     I_1=I_0;
%     if I_dimension==3
%         if image_size(3)==4;
%             I_1=I_1(:,:,1:3);
%             if dynamic_class_0==2^16;
%                 I_1(I_1>=2500)=2500;  % for brightfield image
%                 imagesc(I_1)
%                 colormap(gray)
%             else
%                 imshow(I_1)
%             end
%         elseif image_size(3)==3;
%             if dynamic_class_0==2^16;
%                 I_1(I_1>=2500)=2500;  % for brightfield image
%                 imagesc(I_1)
%                 colormap(gray)
%             else
%                 imshow(I_1)
%             end
%         end
%     else
%         if dynamic_class_0==2^16;
%             I_1(I_1>=2500)=2500;  % for brightfield image
%             imagesc(I_1)
%             colormap(gray)
%         else
%             imshow(I_1)
%         end
%     end
%     axis image
%     axis off
%     hold on
%     plot(thisBoundary_y, thisBoundary_x,'b', 'LineWidth', 3)
%     saveas(gcf,[filesname_1_1 filesname_1 '_with_ROI_using_optiaml_parameters.jpg']);  % 圖片存檔
%     
%     % Performance evaluation of edge edtection
%     cd(codepath);
%     [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
%     cd(filepath);
%     Dice_optimal_parameter(i)=Dice;
%     Jaccard_optimal_parameter(i)=Jaccard;
%     
%     location_optimal_all_mean_1(i)=location_optimal_all_mean;   
    
    dlmwrite(['fitting_points_'  filesname_1_1 '_' num2str(C_sigma^2) '_'  num2str(C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(2)-C_high_threshold_1(1)) '_' num2str(C_high_threshold_1(end)) '_sigmoid_4_parameters.txt'],cf_points_sigmoid_4_parameters);
    dlmwrite(['fitting_points_'  filesname_1_1 '_sigmoid_4_parameters.txt'],cf_points_sigmoid_4_parameters);
    dlmwrite(['fitting_slope_'  filesname_1_1 '_sigmoid_4_parameters.txt'],fx_sigmoid_4_parameters_1);
    dlmwrite(['fitting_slope_slope_' filesname_1_1 '_sigmoid_4_parameters.txt'],fxx_sigmoid_4_parameters_1); 
    dlmwrite(['fitting_r2_'  filesname_1_1 '_sigmoid_4_parameters.txt'],distance_fitting_r2_sigmoid_4_parameters);
    dlmwrite(['locations_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_min_fitting_slope_sigmoid_4_parameters);
    dlmwrite(['locations_mirror_'  filesname_1_1 '_sigmoid_4_parameters.txt'],location_mirror_min_fitting_slope_sigmoid_4_parameters); 
    dlmwrite(['fitting_points_' filesname_1_1 '_' poly_number_1 '.txt'],cf_points_poly);
    dlmwrite(['fitting_slope_'  filesname_1_1 '_'  poly_number_1 '.txt'],fx_poly_1);
    dlmwrite(['fitting_slope_slope_'  filesname_1_1 '_'  poly_number_1 '.txt'],fxx_poly_1); 
    dlmwrite(['fitting_r2_'  filesname_1_1 '_' poly_number_1 '.txt'],distance_fitting_r2_poly); 
    dlmwrite(['segmentation_map_optimal_'  filesname_1_1 '.txt'],bw);
    dlmwrite(['segmentation_edge_optimal_'  filesname_1_1 '.txt'],bw_edge);
    dlmwrite(['locations_'  filesname_1_1 '_'  poly_number_1 '.txt'],location_appropriate_intersection_distance_fitting); 
    estimation_time(i)=toc;
    
%     % Performance evaluation for Matlab default
%     clear thisBoundary thisBoundary_x thisBoundary_y
%     [BW,thres] = edge(I,'Canny');
%     sel1=strel('square',3);
%     sel2=strel('square',5);
%     BW=imdilate(BW,sel2);
%     BW = imfill(BW, 'holes');
%     imLabel = bwlabel(BW);% ??通?域?行?
%     stats = regionprops(imLabel,'Area');
%     [b,index]=sort([stats.Area],'descend');
%     if length(stats)<size_order_threshold  % 前幾大
%         bw=imLabel;
%     else
%         bw=ismember(imLabel,index(1:size_order_threshold));  % 前幾大
%     end
%     cd(codepath);
%     [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores( bw_manual_segmentation, bw);
%     cd(filepath);
%     Dice_Matlab_default(i)=Dice;  
%     Jaccard_Matlab_default(i)=Jaccard; 
%     bw_edge=bwboundaries(bw);
%     k=size(bw_edge,1);
%     thisBoundary = bw_edge{k};
%     thisBoundary_y=thisBoundary(:,2);
%     thisBoundary_x=thisBoundary(:,1);
%     figure('outerposition',get(0,'screensize'));  % 設定圖與螢幕一樣大
%     I_1=I_0;
%     if I_dimension==3
%         if image_size(3)==4;
%             I_1=I_1(:,:,1:3);
%             if dynamic_class_0==2^16;
%                 I_1(I_1>=2500)=2500;  % for brightfield image
%                 imagesc(I_1)
%                 colormap(gray)
%             else
%                 imshow(I_1)
%             end
%         elseif image_size(3)==3;
%             if dynamic_class_0==2^16;
%                 I_1(I_1>=2500)=2500;  % for brightfield image
%                 imagesc(I_1)
%                 colormap(gray)
%             else
%                 imshow(I_1)
%             end
%         end
%     else
%         if dynamic_class_0==2^16;
%             I_1(I_1>=2500)=2500;  % for brightfield image
%             imagesc(I_1)
%             colormap(gray)
%         else
%             imshow(I_1)
%         end
%     end
%     axis image
%     axis off
%     hold on
%     plot(thisBoundary_y, thisBoundary_x,'b', 'LineWidth', 3)
%     saveas(gcf,[filesname_1_1 '_with_ROI_Matlab_default.jpg']);  % 圖片存檔
%     
%     dlmwrite(['segmentation_map_Matlab_default_'  filesname_1_1 '.txt'],bw);
%     dlmwrite(['segmentation_edge_Matlab_default_'  filesname_1_1 '.txt'],bw_edge);
    
    close all
end

dlmwrite(['Estimation_time.txt'],estimation_time);
dlmwrite(['Entropy_of_whole_image_for_threshold_range_selection.txt'],hs_sum_counts);
% dlmwrite(['location_optimal.txt'], location_optimal_all_mean_1);
% dlmwrite(['location_Dice_index_max.txt'],[Dice_index_max(:,2:3) distance_1_max_Dice]);
% dlmwrite(['location_Jaccard_index_max.txt'],[Jaccard_index_max(:,2:3) distance_1_max_Jaccard]);
dlmwrite(['Dice_Jaccard_index_max.txt'],[Dice_index_max(:,1) Jaccard_index_max(:,1)]);
% dlmwrite(['Dice_Jaccard_index_optimal_parameters.txt'],[Dice_optimal_parameter;Jaccard_optimal_parameter]');
% dlmwrite(['Dice_Jaccard_index_Matlab_defaults.txt'],[Dice_Matlab_default;Jaccard_Matlab_default]');
