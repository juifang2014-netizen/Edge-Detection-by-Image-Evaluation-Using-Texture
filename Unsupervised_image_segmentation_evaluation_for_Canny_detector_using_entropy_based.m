clear all
clc
close all
warning off

filepath='C:\Users\User\Desktop\REVISE 1\Testing image\Type_XII';
codepath='C:\Users\User\Desktop\程式';
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

C_threshold_number=2^8;
C_threshold_number_threshold_0=0.95;
outlier_selection=1; % 去除 outluers: 1，不去除 outliers: 2
C_threshold_number_threshold=1-C_threshold_number_threshold_0;
size_order_threshold=1;
lowHighRatio = 0.33;  % 若不給則自適應，低閾值 / 高閾值比例

kk=1;
for i=1:lengthFiles;
    % 每張圖的 threshold range 不同，需要清除相關參數
    tic;
    cd(filepath);
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
    
    dynamic_class_0=class(I_0);
    if strcmp(dynamic_class_0,'uint16')==1
        dynamic_class_0=2^16;
    elseif strcmp(dynamic_class_0,'uint8')==1
        dynamic_class_0=2^8;
    else
    end
    
    % 輸入手動分割的結果
    bw_manual_segmentation=load(['RANGE_' filesname_1_1 '.txt']);
    
     % entropy-based adaptive threshold detection
    cd(codepath);
    [C_low_threshold, C_high_threshold] = adaptiveCannyThreshold_entropy(I);
    cd(filepath);        
    BW = edge(I,'Canny',[C_low_threshold C_high_threshold], C_sigma);   
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
    saveas(gcf,[filesname_1_1 '_entropy_with_ROI.jpg']);  % 圖片存檔
         
    % Performance evaluation of edge edtection
    cd(codepath);
    [Accuracy, Sensitivity, Fmeasure, Precision, MCC, Dice, Jaccard, Specitivity] = EvaluateImageSegmentationScores(bw_manual_segmentation, bw);
    cd(filepath);
    Dice_1(kk)=Dice;
    Jaccard_1(kk)=Jaccard;

    estimation_time(kk)=toc;
    kk=kk+1;
    close all    
end

dlmwrite(['Jaccard_index_entropy.txt'],Jaccard_1);
dlmwrite(['Dice_index_entropy.txt'],Dice_1);
dlmwrite(['Estimation_time_entropy.txt'],estimation_time);
