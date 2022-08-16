% Shell_SF_heterogeneity.m
% This script takes in an image file of fibers and calculates the
% heterogeneity of alignment of the fibers based on methods by:
%   (Fomovsky & Holmes 2010 AJP Heart Circ Physiol)
%   (Richardson & Holmes 2016 Biophys J)
%
% Thien-Khoi N. Phung (August 16, 2022)

%% Calculate Stress Fiber Vectors
    % Load SF file
    SFfile  = 'test_images/U13_72_P1_5_SF.tif';
    SFimage = imread(SFfile);
    SFinfo  = imfinfo(SFfile);
    pr      = 1/SFinfo.XResolution;

    % Load basal cell segmentation
    load('test_images/U13_72_P1_5_basalcell.mat')
    
    % Parameters for Analysis
    % Window side length (um)
    % windowside = 1.5;
    windowside = 3;
    
    % Dilate basal boundaries for exclusion of thick basal borders
    % Create disk with radius of 0.75 um
    disq = strel('disk',ceil(0.75/pr));
    % Erode borders
    basalbodyerode = imerode(basalbody,disq);
    
    % Generate pixel intensity threshold based on basal cell body
    % segmentation superimposed on SF image
    % Using Basal Cell Segmentation
    % nth percentile of pixel data based on cell body
    pthresh =  prctile(SFimage(basalbodyerode),35);
    
    % Border overlap threshold (what % of window pixels can overlap
    % non-basal cell border?
    bthresh = 0;

    % Generate all stress fibers 
    [sfxyai,bcount,sfxyaiall]  = sfvectors(SFimage,...
                                     'window',windowside,...
                                     'pixres',pr,...
                                     'threshold',pthresh,...
                                     'borderoverlap',~basalbodyerode,...
                                     'borderthresh',bthresh,...
                                     'visualize',true);
    
%% Stress Fiber Heterogeneity Metric
    SFXYAIB = [sfxyaiall bcount];
    % Which of the SF to keep in the analysis
    % No basal cell border overlap in the window & intensity>pthresh
    keepvecs             = SFXYAIB(:,5)==0 & SFXYAIB(:,4)>=pthresh;
    SFXYAIB(~keepvecs,:) = [];
        
    [hetauc,hetmet]   = sfhetmetric(SFXYAIB,pr);

    figure('NumberTitle','off','Name','HetMet',...
           'Units','Inches','Position',[3 3 5 2.5])
    ax = axes('NextPlot','add',...
              'FontName','Open Sans','FontSize',10,...
              'FontWeight','bold',...
              'TickDir','out',...
              'LineWidth',1.5,...
              'Clipping','on');
    ff = fill([hetmet(:,1); flipud(hetmet(:,1))],...
         [hetmet(:,2);flipud(hetmet(:,4))],'r');
    ff.FaceColor = [239,138,98]./255;
    ff.EdgeColor = 'none';    
    hold on
    plot(hetmet(:,1),hetmet(:,2),'r-','LineWidth',3)
    plot(hetmet(:,1),hetmet(:,4),'r:','LineWidth',3)
    
    xlim([0 100])
    xticks([0 100])
    xlabel('Distance Between Fibers')
    ylim([0.45 0.75])
    yticks([0.45 0.75])
    ylabel('Mean Fiber-Fiber Alignment')

    legend('AUC','Experimental Data','Randomized (control)')