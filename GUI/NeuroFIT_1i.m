function [varargout] = NeuroFIT_1i(GUIfh,varargin)
%% NeuroFIT - Neuro FLUORESCENT IMAGING TOOLBOX



%% -- SETTING UP MAIN FIGURE & MESSAGE CON
clc; close all; %clear all;

Fh1 = figure(1);
set(Fh1,'OuterPosition',[550 400 1100 700],'Color',[1,1,1],'Tag','GUIfh')
hax1 = axes('Position',[.05 .22 .44 .7],'Color','none','XTick',[],'YTick',[]);
hax2 = axes('Position',[.53 .22 .44 .7],'Color','none','XTick',[],'YTick',[]);
% hax3 = axes('Position',[.05 .02 .9 .15],'Color','none','XTick',[],'YTick',[]);
% axes(hax1); set(Fh1,'CurrentAxes',hax1); get(hax3); %clf(GUIfh)


spfN=sprintf(' ');spf1=sprintf('>>'); spf2=sprintf('>>');spf3=sprintf('>>');spf4=sprintf('>>');
str = {spfN, spf1,spf2,spf3,spf4};
ft = annotation(Fh1,'textbox', [0.05,0.03,0.7,0.13],'String', str,'FontSize',14);
set(ft,'interpreter','none') 


%% -- DEAL INPUT ARGS

    if exist('varargin','var') && nargin > 0
        [STRs, DODs, NUMs, GHAXs, DOTs] = deal(varargin{:});
    else
        %--- PRINT MESSAGE TO CON ---
        spf0=sprintf('Aborting non-GUI run attempt');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        spf0=sprintf('Use NeuroFIT_GUI.m to run this function');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        spf0=sprintf('(or use the GUI that just opened)');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        pause(3)
        %----------------------------
        return
    end



%% -- GET DECISION TREE INFO
%---------------------------------------------

ImageFile1 = STRs{1};
ImageFile2 = STRs{2};
DataSaveFileName = STRs{3};
ZstackSaveFileName = STRs{4};

doPromptMediaDirName1 = DOTs(12);
doPromptMediaDirName2 = DOTs(13);

dosavedata = DOTs(14);
doCSV = DOTs(15);
doMAT = DOTs(16);


%% -- GET ALL THE THE MEDIA FILES LOCATED IN SPECIFIED FOLDER

if doPromptMediaDirName1

    [iFileName,iPathName] = uigetfile('*.tif*','Select image file');

else

iFileName = ImageFile1;

%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Getting image file(s): %s',ImageFile1);
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4); pause(.2)
%----------------------------
end




%% -- RESIZE IMAGES AND STORE IMAGE PIXEL MATRIX VALUES INTO A CELL ARRAY

Pixels = [512 NaN];         % resize all images to 512x512 pixels

%--- PRINT MESSAGE TO CON ---
spf0=sprintf('resizing file(s) to: % 4.4g x % 4.4g pixels',Pixels(1),Pixels(1));
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    [I,map] = imread(iFileName);   % get image data from file

    % colormap will be a 512x512 matrix of class double (values range from 0-1)
    iDUBs = im2double(I);              
    iDUBs = imresize(iDUBs, Pixels);
    iDUBs(iDUBs > 1) = 1;  % In rare cases resizing results in some pixel vals > 1
    
    szImg = size(iDUBs);

if numel(szImg) == 3
    AveR = mean(mean(iDUBs(:,:,1)));
    AveG = mean(mean(iDUBs(:,:,2)));
    AveB = mean(mean(iDUBs(:,:,3)));

    if AveR > AveG
        iDUB = iDUBs(:,:,1);
        %--- PRINT MESSAGE TO CON ---
        spf0=sprintf('pixel intensities strongest in red channel');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        %----------------------------
    else
        iDUB = iDUBs(:,:,2);
        %--- PRINT MESSAGE TO CON ---
        spf0=sprintf('pixel intensities strongest in green channel');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        %----------------------------
    end

        axes(hax1);
    imagesc(iDUBs(:,:,1)); title('red channel');
        axes(hax2);
    imagesc(iDUBs(:,:,2)); title('green channel');

end

pause(1)

%% -- NORMALIZE DATA TO RANGE: [0 <= DATA <= 1]
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Normalizing data to range [0 <= DATA <= 1]');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------


IMGsumOrig = iDUB;
IMGs = iDUB;

maxIMG = max(max(IMGs));
minIMG = min(min(IMGs));

lintrans = @(x,a,b,c,d) (c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)));

    for nn = 1:numel(IMGs)

        x = IMGs(nn);

        if maxIMG > 1
            IMGs(nn) = lintrans(x,minIMG,maxIMG,0,1);
        else
            IMGs(nn) = lintrans(x,minIMG,1,0,1);
        end


    end


%% ----- VIEW ORIGINAL VS NORMALIZED IMAGES  -----
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Displaying original vs normalized image');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    axes(hax1); %axes(GUIfh.Children(1));
imagesc(IMGsumOrig); title('Pre-Normalized Image')

    axes(hax2); %axes(Fh1.Children(2));
imagesc(IMGs); title('Normalized Image')
    
    colormap('bone')

pause(2)






%% USE MOUSE TO DRAW BOX AROUND BACKGROUND AREA

if DODs(2)
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Perform manual background selection...');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    iDUB = IMGs;

        fh11 = figure(11);
        set(fh11,'OuterPosition',[400 400 700 700])
        ax1 = axes('Position',[.1 .1 .8 .8]);
    imagesc(iDUB);
        title('USE MOUSE TO DRAW BOX AROUND BACKGROUND AREA')
        % colormap(bone)

        disp('DRAW BOX AROUND A BACKGROUND AREA')
    h1 = imrect;
    pos1 = round(getPosition(h1)); % [xmin ymin width height]

close(fh11)

figure(Fh1)
set(Fh1,'CurrentAxes',hax1);
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Perform manual background selection... done');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
end


if DODs(3)
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Performing auto background selection');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    iDUB = IMGs;
    szBG = size(iDUB);
    BGrows = szBG(1);
    BGcols = szBG(2);
    BGr10 = floor(BGrows/10);
    BGc10 = floor(BGcols/10);
    pos1 = [BGrows-BGr10 BGcols-BGc10 BGr10-1 BGc10-1];

end



    
%% GET FRAME COORDINATES AND CREATE XY MASK

    pmsk = [pos1(2) (pos1(2)+pos1(4)) pos1(1) (pos1(1)+pos1(3))];

    % Background
    mask{1} = zeros(size(iDUB));
    mask{1}(pmsk(1):pmsk(2), pmsk(3):pmsk(4)) = 1;
    mask1 = mask{1};


%% CHECK THAT MASK(S) ARE CORRECT

        axes(hax1);
    imagesc(iDUB);
        
vert = [pmsk(1) pmsk(1);pmsk(1) pmsk(2);pmsk(2) pmsk(2);pmsk(2) pmsk(1)];
fac = [1 2 3 4]; % vertices to connect to make square
fvc = [1 0 0;0 1 0;0 0 1;0 0 0];
patch('Faces',fac,'Vertices',vert,...
'FaceVertexCData',fvc,'FaceColor','interp','FaceAlpha',.5)

pause(1)


%% -- GET MEAN OF BACKGROUND PIXELS & SUBTRACT FROM IMAGE
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Taking average of background pixels');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    f1BACKGROUND = iDUB .* mask1;
    meanBG = mean(f1BACKGROUND(f1BACKGROUND >= 0));

    meanALL = mean(iDUB(:));

    iDUB = iDUB - meanBG;
    iDUB(iDUB <= 0) = 0;



%% -- MESH SURFACE PLOT
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('displaying image intensity peaks as 3D mesh');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

        axes(hax2);
    mesh(iDUB)
        view([186 24])
        title('image intensity peaks as 3D mesh')
        pause(2)

%% check vague SNR for masks.
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('plotting histogram of pixel intensities');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    hist(hax1,iDUB(:),80);
        title('histogram of pixel intensities')
        pause(.5)


%% SET THRESHOLD BASED ON SNR HISTOGRAM

if DODs(4)
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('manually set threshold for target inclusion...');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    promptTXT = {'Enter Threshold Mask Values:'};
    dlg_title = 'Input'; num_lines = 1; 

    presetval = {num2str(NUMs(1))};
    
    dlgOut = inputdlg(promptTXT,dlg_title,num_lines,presetval);
    threshmask = str2num(dlgOut{:});

%--- PRINT MESSAGE TO CON ---
spf0=sprintf('manually set threshold for target inclusion...done');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
end


if DODs(5)
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('using preset threshold for target inclusion');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    threshmask = NUMs(1);

end


%% -- REMOVE PIXELS BELOW THRESHOLD
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('removing pixels below threshold');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    threshPix = iDUB > threshmask;  % logical Mx of pixels > thresh
    rawPix = iDUB .* threshPix;		% raw value Mx of pixels > thresh


%% -- CHOOSE PCT% RANGE OF PIXELS ABOVE THRESH

    % how many pixels passed threshold (tons!)?  
	n = sum(threshPix(:));

    % get actual values of those pixels
	valArray = iDUB(threshPix);

    % sort pixels, brightest to dimmest
	Hi2LoVals = sort(valArray, 'descend');


%% -- CHOOSE PCT% RANGE OF PIXELS ABOVE THRESH

if DODs(4)
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('waiting for entry of target pixels analysis range...');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

	promptTxtUB = {'Enter upper-bound percent of pixels to analyze'};
	dlg_TitleUB = 'Input'; num_lines = 1; presetUBval = {'99.99'};
	UB = inputdlg(promptTxtUB,dlg_TitleUB,num_lines,presetUBval);
	UpperBound = str2double(UB{:}) / 100;

	promptTxtLB = {'Enter lower-bound percent of pixels to analyze'};
	dlg_TitleLB = 'Input'; num_lines = 1; presetLBval = {'5'};
	LB = inputdlg(promptTxtLB,dlg_TitleLB,num_lines,presetLBval);
	LowerBound = str2double(LB{:}) / 100;


%--- PRINT MESSAGE TO CON ---
spf0=sprintf('waiting for entry of target pixels analysis range...done');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
else
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('using preset range for target pixel analysis');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    UpperBound = NUMs(2)/ 100;
    LowerBound = NUMs(3)/ 100;

end

	n90 = round(n - (n * UpperBound));
    if n90 < 1; n90=1; end;
	n80 = round(n - (n * LowerBound));
	hotpix = Hi2LoVals(n90:n80);
	

%% -- GET PIXELS THAT PASSED PCT% THRESHOLD
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('getting pixels inside analysis range');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    HighestP = Hi2LoVals(n90);
    LowestsP = Hi2LoVals(n80);

    HiLogicMxP = iDUB <= HighestP;			% logic value Mx of pixels passed thresh
    HiRawMxP = iDUB .* HiLogicMxP;		% raw value Mx of pixels passed thresh
    LoLogicMxP = iDUB >= LowestsP;			% logic value Mx of pixels passed thresh
    LoRawMxP = iDUB .* LoLogicMxP;		% raw value Mx of pixels passed thresh

    IncLogicMxP = HiRawMxP > LowestsP;
    IncRawMxP = HiRawMxP .* IncLogicMxP;


    IncPixArray = IncRawMxP(IncRawMxP>0);
    Hi2LoIncPixArray = sort(IncPixArray, 'descend');


%% -- GET MEAN MEDIAN & SD OF REMAINING PIXELS
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('computing stats for pixels inside analysis range');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
	
	mt = mean(Hi2LoIncPixArray);
	mdn = median(Hi2LoIncPixArray);
	st = std(Hi2LoIncPixArray);

%--- PRINT MESSAGE TO CON ---
spf0=sprintf('Mean: % 4.4g ',mt);
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
spf0=sprintf('Median: % 4.4g ',mdn);
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
spf0=sprintf('StDev: % 4.4g ',st);
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
pause(2)
%----------------------------


%% -- Boxplot & Histogram
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('creating boxplot and histogram of analyzed pixels');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    boxplot(hax1,Hi2LoIncPixArray ...
	,'notch','on' ...
	,'whisker',1 ...
	,'widths',.8 ...
	,'factorgap',[0] ...
	,'medianstyle','target');
        title('boxplot of pixels included in analysis')
        pause(1)


    hist(hax2,Hi2LoIncPixArray(:),100);
        title('histogram of pixels included in analysis')
		pause(3)


%% -- VIEW ORIGINAL, SNR MASK, AND TARGET IMAGES
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('showing original, snr mask, and target images');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    %-------
    %fig23 = figure(23);
    %set(fig23,'OuterPosition',[10  50  1600  600],'Color',[1,1,1])
    %-------

        %HaxOrig = axes('Position',   [.03 .05 .3 .9]);
        axes(hax1);
    imagesc(IMGsumOrig);
        title('Original Image')

        %HaxThresh = axes('Position', [.36 .05 .3 .9]);
    %imagesc(iDUB);

        %HaxInc = axes('Position',    [.69 .05 .3 .9]);
        axes(hax2);
    imagesc(IncRawMxP);
         title('Pixels targeted for analysis')

    colormap('bone');

%% --- SET COLORMAP --

    cmaps = cmaplist;
    % [smap,butn] = listdlg('PromptString','Select a Colormap:',...
    % 'SelectionMode','single','ListSize',[200 200],'ListString',cmaps);

   smap = 1;
    if smap==1; cmap=customcmap(1); else cmap=cmaps{smap}; end
    colormap(cmap); pause(2); %close gcf;

    gmap=customcmap(2);
    rmap=customcmap(1);


%% --- PERFORM CONVOLUTION WITH GAUSSIAN FILTER ---
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('performing pixel smoothing via gaussian matrix convolution');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    doConvnFilter = 1;
    if doConvnFilter

        % MaskMaker([PEAK HEIGHT] [MASK SIZE] [SLOPE SD] [RESOLUTION])
        Mask = MaskMaker(2.5, 5, .17, .1);
        fIMG = convn(IncRawMxP,Mask,'same');

        %axes(Fh1.Children(1).Children(2))
        set(Fh1,'CurrentAxes',hax1);
        imagesc(fIMG)
            title('smoothed image')
        colormap(cmap)
            

    else
        fIMG = IncRawMxP;
    end


%% --- PERFORM CELL AUTOTRACE ---
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('performing cell autotracing');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4); pause(.2)
%----------------------------

    % Pad the image with zeros so no cell "edge" touches the image outer boarder
    fIMG(1:5,:) = 0;
    fIMG(end-5:end,:) = 0;
    fIMG(:,1:5) = 0;
    fIMG(:,end-5:end) = 0;

    dim = size(fIMG);
    col = round(dim(2)/2)-90;
    row = min(find(fIMG(:,col)));

    fIMG_filled = imfill(fIMG,'holes');
    fbIMG = bwboundaries(fIMG_filled);

    for k=1:numel(fbIMG)
       fbkIMG(k) = numel(fbIMG{k});
    end
    

fbkThresh = 60;
    tIMG = fbkIMG < fbkThresh;
    fbIMG(tIMG) = [];



%% -- ASSIGN AXES HANDLES

% cla(Fh1.Children(1).Children(2))
% cla(Fh1.Children(1).Children(1))
% hax1 = Fh1.Children(1).Children(1);
% hax2 = Fh1.Children(1).Children(2);


%% --- PLOT AUTOTRACED EDGES ---
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('plotting cell autotraces');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
cla(hax1)
cla(hax2)

        axes(hax1)
    imagesc(fIMG)
        colormap(hax1,gmap)
        title('Convolution Processed Image')
        hold on
        

        for k=1:numel(fbIMG)
            plot(hax1,fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
        end
        hold on

pause(2)

        axes(hax2)
    imagesc(IMGsumOrig)
        colormap(hax2,gmap)
        title('Original Image Sum')
        hold on

    for k=1:numel(fbIMG)
        plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w','LineWidth',2);
    end
        hold on

pause(2)


%% --- GET IMAGE TO COMPARE AUTOTRACED BETWEEN CHANNELS ---
    
    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('getting complementary color channel image');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

    iDUBg = IMGsumOrig;

% promptTXT = {'Enter Image Filename'};
%     dlg_title = 'Input'; num_lines = 1; presetval = {'NFIT_IMG_512px_RGB_16bit.tif'};
%     dlgOut = inputdlg(promptTXT,dlg_title,num_lines,presetval);
%     imgname = dlgOut{:};


    rFileName = 'FLred.tif';

    [f2,map] = imread(rFileName);		% import image
    iDUBsr = im2double(f2);             % colormap is Mx-by-1 array of class double.
    iDUBsr = imresize(iDUBsr, [512 NaN]);
    iDUBsr(iDUBsr > 1) = 1;


%% -- CONVERT SECOND IMAGE TO USABLE MATRIX
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('processing complementary image');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
cla(hax1)
cla(hax2)

    szImg = size(iDUBsr);

if numel(szImg) == 3
    AveR = mean(mean(iDUBsr(:,:,1)));
    AveG = mean(mean(iDUBsr(:,:,2)));
    AveB = mean(mean(iDUBsr(:,:,3)));

    if AveR > AveG
        iDUBr = iDUBsr(:,:,1);
        %--- PRINT MESSAGE TO CON ---
        spf0=sprintf('pixel intensities strongest in red channel');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        %----------------------------
    else
        iDUBr = iDUBsr(:,:,2);
        %--- PRINT MESSAGE TO CON ---
        spf0=sprintf('pixel intensities strongest in green channel');
        [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
        %----------------------------
    end

end

    szIMG = size(iDUBr);

    if numel(szIMG) > 2
        if szIMG(3) == 3
                Im = rgb2gray(iDUBr);
        elseif szIMG(3) == 4
                Im = iDUBr(:,:,1);
        else
                Im = iDUBr(:,:,1);
        end
      iDUBr = Im;
    end

        axes(hax2)
    imagesc(iDUBr)
        title('Comparison Image')
        colormap(hax2,rmap)


%% -- DISPLAY COMPLIMENTARY IMAGE
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('displaying complimantary images');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
cla(hax1)
cla(hax2)

        axes(hax1)
    imagesc(iDUBg)
        colormap(hax1,gmap)
        title('Green Channel Processed Image')
        hold on
        

        for k=1:numel(fbIMG)
    plot(hax1,fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
        end
        hold on


pause(.5)

        axes(hax2)
    imagesc(iDUBr)
        colormap(hax2,rmap)
        title('Red Channel Comparison Image')
        hold on

    for k=1:numel(fbIMG)
        plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w','LineWidth',2);
    end
        hold on

pause(2)



%% -- CREATE BINARY IMAGE
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('creating binary image');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

    % Create black and white (BW) binary image
    fIMG_BW = im2bw(fIMG_filled, .5);


%% -- ASSIGN LABELS TO OBJECTS IN BINARY IMAGE
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('assigning labels to discrete objects in binary image');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
cla(hax1)
cla(hax2)

    % Returns labels for connected components of a binary image
    [L, NUM] = bwlabeln(fIMG_BW);  

        axes(hax1)
    imagesc(fIMG_BW);
        title('Binary Image')

        axes(hax2)
    imagesc(L);
        title('Binary Objects')

pause(3)
        

%% -- CALCULATE STATISTICS FROM REGIONPROPS AT LOCATIONS (L) ON IncRawMxP

%--- PRINT MESSAGE TO CON ---
spf0=sprintf('calculating pixel intensity statistics for target objects');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

STATS = regionprops(L,IncRawMxP, 'MeanIntensity','Perimeter','PixelList','PixelValues')
%STATS = regionprops(L,IncRawMxP, 'all')



    for n = 1:size(STATS,1)
        perim(n) = STATS(n).Perimeter >= (fbkThresh/2);
    end

    for n = 1:size(STATS,1)
        PixIntensity(n) = STATS(n).MeanIntensity * perim(n);
    end

    [PIrow,PIcol,PixelVals] = find(PixIntensity)
    nPixIntensity = numel(PIrow);


    for n = 1:size(STATS,1)
        PixList{n} = STATS(n).PixelList;
    end
    PixList(~perim) = [];



%% -- CHECK THAT STATS.PixelList CONTAINS CORRECT PIXEL LOCATIONS
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('determining if stats were crunched for correct pixels');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------

PLBW = zeros(size(fIMG_BW));

for nn = 1:numel(PixList)

    PLnn = PixList{nn}';

    for mm = 1:numel(PLnn(1,:))

        PLBW(PLnn(2,mm),PLnn(1,mm)) = 1;

    end

end



%% -- DO SOME NORMALIZATION ON STATS.PixelValues
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('normalizing data');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------


PVs = PixelVals;

maxPV = max(max(PVs));
minPV = min(min(PVs));

if maxPV > 1
lintrans = @(x,a,b,c,d) (c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)));

    for nn = 1:numel(PVs)

        x = PVs(nn);

        PVs(nn) = lintrans(x,minPV,maxPV,0,1);

    end

end

%% -- PLOT BAR GRAPH OF PIXEL INTENSITIES FOR TRACED CELLS
%--- PRINT MESSAGE TO CON ---
spf0=sprintf('plotting bar graph of pixel intensities for traced cells');
[ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------
cla(hax1)
cla(hax2)
	
    axes(hax1), hold off
bar(PVs)
    title('Pixel intensities for traced cells')
	
    axes(hax2), hold off
imagesc(PLBW);
    title('Traced cells object locations')


pause(1)

finalIMG = iDUB;
finalIMG(~PLBW) = 0;

        axes(hax2)
    imagesc(finalIMG)
        colormap(hax2,gmap)
        title('Final Cell Trace')
        hold on
        

        for k=1:numel(fbIMG)
    plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w:','LineWidth',1.5);
        end
        hold on

pause(2)



%% -- END FUNCTION
varargout = {mt,mdn,st};
end
%====================================================================






function ccmap = customcmap(n)

switch n

    case 1

        % mCherry colormap
        ccmap = [
        0.0417         0         0;
        0.0833         0         0;
        0.1250         0         0;
        0.1667         0         0;
        0.2083         0         0;
        0.2500         0         0;
        0.2917         0         0;
        0.3333         0         0;
        0.3750         0         0;
        0.4167         0         0;
        0.4583         0         0;
        0.5000         0         0;
        0.5417         0         0;
        0.5833         0         0;
        0.6250         0         0;
        0.6667         0         0;
        0.7083         0         0;
        0.7500         0         0;
        0.7917         0         0;
        0.8333         0         0;
        0.8750         0         0;
        0.9167         0         0;
        0.9583         0         0;
        1.0000         0         0;
        1.0000    0.0417         0;
        1.0000    0.0833         0;
        1.0000    0.1250         0;
        1.0000    0.1667         0;
        1.0000    0.2083         0;
        1.0000    0.2500         0;
        1.0000    0.2917         0;
        1.0000    0.3333         0;
        1.0000    0.3750         0;
        1.0000    0.4167         0;
        1.0000    0.4583         0;
        1.0000    0.5000         0;
        1.0000    0.5417         0;
        1.0000    0.5833         0;
        1.0000    0.6250         0;
        1.0000    0.6667         0;
        1.0000    0.7083         0;
        1.0000    0.7500         0;
        1.0000    0.7917         0;
        1.0000    0.8333         0;
        1.0000    0.8750         0;
        1.0000    0.9167         0;
        1.0000    0.9583         0;
        1.0000    1.0000         0;
        1.0000    1.0000    0.0625;
        1.0000    1.0000    0.1250;
        1.0000    1.0000    0.1875;
        1.0000    1.0000    0.2500;
        1.0000    1.0000    0.3125;
        1.0000    1.0000    0.3750;
        1.0000    1.0000    0.4375;
        1.0000    1.0000    0.5000;
        1.0000    1.0000    0.5625;
        1.0000    1.0000    0.6250;
        1.0000    1.0000    0.6875;
        1.0000    1.0000    0.7500;
        1.0000    1.0000    0.8125;
        1.0000    1.0000    0.8750;
        1.0000    1.0000    0.9375;
        1.0000    1.0000    1.0000];

    case 2

    % GFP colormap
    ccmap =[
         0         0         0;
         0    0.0389         0;
         0    0.0778         0;
         0    0.1167         0;
         0    0.1556         0;
         0    0.1945         0;
         0    0.2335         0;
         0    0.2685         0;
         0    0.3035         0;
         0    0.3385         0;
         0    0.3735         0;
         0    0.4150         0;
         0    0.4565         0;
         0    0.4980         0;
         0    0.5367         0;
         0    0.5753         0;
         0    0.6139         0;
         0    0.6525         0;
         0    0.6911         0;
         0    0.7297         0;
         0    0.7683         0;
         0    0.8069         0;
         0    0.8456         0;
         0    0.8842         0;
         0    0.9228         0;
         0    0.9614         0;
         0    1.0000         0;
    0.0286    0.9917    0.0087;
    0.0573    0.9834    0.0175;
    0.0859    0.9751    0.0262;
    0.1145    0.9668    0.0350;
    0.1432    0.9585    0.0437;
    0.1718    0.9502    0.0525;
    0.2005    0.9419    0.0612;
    0.2291    0.9336    0.0700;
    0.2577    0.9253    0.0787;
    0.2864    0.9170    0.0875;
    0.3150    0.9087    0.0962;
    0.3436    0.9004    0.1049;
    0.3723    0.8921    0.1137;
    0.4009    0.8838    0.1224;
    0.4295    0.8755    0.1312;
    0.4826    0.8801    0.1275;
    0.5356    0.8848    0.1239;
    0.5886    0.8894    0.1203;
    0.6416    0.8940    0.1166;
    0.6947    0.8986    0.1130;
    0.7477    0.9033    0.1094;
    0.8007    0.9079    0.1057;
    0.8152    0.9092    0.1047;
    0.8296    0.9104    0.1037;
    0.8441    0.9117    0.1027;
    0.8586    0.9130    0.1017;
    0.8730    0.9142    0.1008;
    0.8875    0.9155    0.0998;
    0.9020    0.9167    0.0988;
    0.9164    0.9180    0.0978;
    0.9309    0.9193    0.0968;
    0.9453    0.9205    0.0958;
    0.9598    0.9218    0.0948;
    0.9699    0.9413    0.2711;
    0.9799    0.9609    0.4474;
    0.9900    0.9805    0.6237;
    1.0000    1.0000    0.8000];

    otherwise

end

end



function cmaps = cmaplist
% -- store colormaps for later use
cmaps = ['custom   ';
         'hot      ';
		 'jet      ';
		 'bone     ';
		 'autumn   ';
		 'colorcube';
		 'cool     ';
		 'copper   ';
		 'flag     ';
		 'gray     ';
		 'hsv      ';
		 'lines    ';
		 'pink     ';
		 'prism    ';
		 'spring   ';
		 'summer   ';
		 'white    ';
		 'winter   '];
cmaps = cellstr(cmaps);
end



function [ft,spf2,spf3,spf4] = upcon(ft,spf0,spf2,spf3,spf4)

spfN=sprintf(' ');
disp(spf0)
spf1=spf2;spf2=spf3;spf3=spf4;spf4=spf0;
ft.String={spfN,spf1,spf2,spf3,spf4};
drawnow; pause(.1)

end




%---------------------------------------------
% NOTES ON HOW THIS WORKS (using MediaDir.m)
%{

Below, we will execute this line of code:

    >> ImageFiles = MediaDir_XXX();

Note:   "_XXX"   will be unique to a specific folder (explained below)


Prior to running the NeuroFIT.m function-script (this file), I created a 
Matlab function called 'MediaDir.m' and put it in the folder containing 
all the desired .tif images we want to analyze during this run (typically 
this will consist of a z-stack of images from a single channel).

I then renamed the file MediaDir.m file to something unique, to identify
the particular folder containing the images. For example, if the folder
containing the images we currently want to analyze is named: "30_35-001"
it would be best practice to make a copy of MediaDir.m and paste it into 
the folder "30_35-001" and then rename the file MediaDir_30_35_001.m
(see 'footNote1' below; also make sure the folder containing the 
images + MediaDir_30_35_001.m are added to your Matlab path).


What does this do?
Not much. When MediaDir_30_35_001() is envoked, matlab will look for a file
called MediaDir_30_35_001.m somewhere in your active path. It will then
evoke the function contained in MediaDir_30_35_001.m -- this function
simply returns all the .tif files located in the same folder as 
MediaDir_30_35_001.m



Is there an easier way to do this?
YES and NO. All we are trying to do is save all the images file names into
a cell array, with one file name per cell. Another way to do this would 
be to navigate to the folder containing the images you want to analyze, 
and run this segment of code:

>> ImageFiles = ls('*.tif*')

The problem with this, is that you now have one long charactar array of all
your file names. Perhaps there is a better way to do all this, but since
the method above works, and there is a lot of other stuff to worry about
I'll leave that for another time.



footNote1: Notice that the folder named "30_35-001" has a dash "-" in its name.
It is bad practice to use dashes or other special charactars in the names
of files and folders. If a special charactar is needed, use the underscore
"_" character. The reason this is bad practice is because a Matlab file:
"MediaDir_30_35-001.m" would be an invalid file name. However, the Matlab
file: "MediaDir_30_35_001.m" would be fine. Generally, 'camel case' is
the preferred way to name files and folders (e.g. camelCaseExample3035001).

%}
%---------------------------------------------
