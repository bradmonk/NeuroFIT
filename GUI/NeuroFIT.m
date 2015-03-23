function [varargout] = NeuroFIT(varargin)
%% NeuroFIT - Neuro FLUORESCENT IMAGING TOOLBOX
clc; %close all; %clear all;

% Change the current folder to the folder of this m-file.
if(~isdeployed); cd(fileparts(which(mfilename))); end


%% -- SETTING UP MAIN FIGURE
% close all
% GUIfh = figure(1);
% set(GUIfh,'OuterPosition',[550 400 1100 800],'Color',[1,1,1],'Tag','GUIfh')
% hax1 = axes('Position',[.05 .22 .42 .7],'Color','none','XTick',[],'YTick',[]);          
% hax2 = axes('Position',[.52 .22 .42 .7],'Color','none','XTick',[],'YTick',[]);
% hax3 = axes('Position',[.05 .02 .9 .15],'Color','none','XTick',[],'YTick',[]);
    % axes(hax1)
    % set(GUIfh,'CurrentAxes',hax1);
    % get(hax3)

%% --- SETTING UP CONSOLE & GUI MESSAGING

GUIfh = figure(NeuroFIT_GUI);figure(GUIfh);

spfN=sprintf(' ');spf1=sprintf('>>'); spf2=sprintf('>>');spf3=sprintf('>>');spf4=sprintf('>>');
str = {spfN, spf1,spf2,spf3,spf4};
ft = annotation(GUIfh,'textbox', [0.05,0.02,0.7,0.1],'String', str,'FontSize',14);
set(ft,'interpreter','none') 

    for x = 1:5;
        spf0 = sprintf('>> launching in... % 6.4g ',6-x);
            [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    end
    for x = 1:5;
        spf0 = sprintf(' ');
            [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    end


%--- PRINT MESSAGE TO CON ---
% spf0=sprintf('TEXTTEXT: % 4.4g ',NUMBER);
% [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
% 
% spf0=sprintf('TEXTTEXTTEXTTEXT');
% [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
%----------------------------


%% --ACCESS COMMAND WINDOW TEXT PROGRAMMATICALLY
% % Get instances of command window text
% jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
% jCmdWin = jDesktop.getClient('Command Window');
% jTextArea = jCmdWin.getComponent(0).getViewport.getView;
% cwText = char(jTextArea.getText);
% % Callback to update GUI whenever text is added to the CW:
% set(jTextArea,'CaretUpdateCallback',@myUpdateFcn)




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


doProcess1imgs = DOTs(1);
doProcess2imgs = DOTs(2);

doProcess1i = DOTs(3);
doProcess1z = DOTs(4);
doProcess1z3D = DOTs(5);


doProcess2iSiDc = DOTs(6);
doProcess2iDiSc = DOTs(7);
doProcess2iDiDc = DOTs(8);
doProcess2zSzDc = DOTs(9);
doProcess2zDzSc = DOTs(10);
doProcess2zDzDc = DOTs(11);

doPromptMediaDirName1 = DOTs(12);
doPromptMediaDirName2 = DOTs(13);
dosavedata = DOTs(14);
doCSV = DOTs(15);
doMAT = DOTs(16);





%% -- REDIRECT TO APPROPRIATE PROCESSING FUNCTION

stats = [];

if doProcess1imgs
    if doProcess1i
        stats = NeuroFIT_1i(GUIfh,STRs,DODs,NUMs,GHAXs,DOTs);
    end
    if doProcess1z
        stats = NeuroFIT_1z(STRs,DODs,NUMs,GHAXs,DOTs)
    end
end

if doProcess2imgs
    if doProcess2iSiDc
        stats = NeuroFIT_2iSiDc(STRs,DODs,NUMs,GHAXs,DOTs);
    end
    if doProcess2iDiSc
        stats = NeuroFIT_2iDiSc(STRs,DODs,NUMs,GHAXs,DOTs)
    end
    if doProcess2iDiDc
        stats = NeuroFIT_2iDiDc(STRs,DODs,NUMs,GHAXs,DOTs)
    end
    if doProcess2zSzDc
        stats = NeuroFIT_2zSzDc(STRs,DODs,NUMs,GHAXs,DOTs)
    end
    if doProcess2zDzSc
        stats = NeuroFIT_2zDzSc(STRs,DODs,NUMs,GHAXs,DOTs)
    end
    if doProcess2zDzDc
        stats = NeuroFIT_2zDzDc(STRs,DODs,NUMs,GHAXs,DOTs)
    end
end

varargout = {stats};
return

try
   return
catch exception
    %rethrow(exception)
end 














%% -- GET ALL THE THE MEDIA FILES LOCATED IN SPECIFIED FOLDER


evStr = strcat('ImageFiles=', STRs{1} ,';');
eval(evStr);    % eval( 'ImageFiles = MediaDir_ExampIMGs_Ch1()' );

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Getting image files using: %s',evStr);
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4); pause(.2)
    %----------------------------






%% -- RESIZE IMAGES AND STORE IMAGE PIXEL MATRIX VALUES INTO A CELL ARRAY

% THE FINAL CELL ARRAY OF INTEREST WILL BE CALLED 'iDUBs'
% To access the pixel value matrix from a single image from this cell array, 
% for example the first image of the z-stack, simply use: 
% >> Img = iDUBs{1}


numF = numel(ImageFiles);   % number of image files
Pixels = [512 NaN];         % resize all images to 512x512 pixels

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('resizing files to: % 4.4g x % 4.4g pixels',Pixels(1),Pixels(1));
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

for nT = 1:numF

    [I,map] = imread(ImageFiles{nT});   % get image data from file
    %infoto = imfinfo(ImageFiles{nT});  % image meta info

    % colormap will be a 512x512 matrix of class uint8 (values range from 0-255)
    iIND = gray2ind(I);               
    iIND = imresize(iIND, Pixels);
    iINDs{nT} = iIND;

    % colormap will be a 512x512 matrix of class double (values range from 0-1)
    iDUB = im2double(I);              
    iDUB = imresize(iDUB, Pixels);
    iDUB(iDUB > 1) = 1;  % In rare cases resizing results in some pixel vals > 1
    iDUBs{nT} = iDUB;

end



%% -- RECORD ANIMATION OF IMAGE Z-STACK (SAVES .avi FILE); REPLAY ANIMATION

if DODs(1)

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Creating Z-stack animation');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

    writerObj = VideoWriter(STRs{4});
    writerObj.FrameRate = 10;
    open(writerObj);


    figure(1)
    set(gcf,'Renderer','zbuffer');
    imagesc(iDUBs{1})
    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';

    for nT = 1:numF

        imagesc(iDUBs{nT})
        title(ImageFiles{nT},'Interpreter','none')
        drawnow

        MovieFrames(nT) = getframe;
        writeVideo(writerObj,MovieFrames(nT));

        pause(.02)
        disp(strcat('iDUBs{',num2str(nT),'}  :: ' ,ImageFiles{nT}))
        pause(.02)

    end

    close(writerObj);
    % figure(2)
    % movie(MovieFrames,2) % To play movie 2x
else

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Skipping Z-stack animation');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

end


%% -- SUM OVER Z-STACK OF IMAGES


    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Flattening z-stack');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

IMGsum = iDUBs{1};

for nT = 2:numF

    IMGx = iDUBs{nT};

    IMGsum = IMGsum + IMGx;

end



%% -- SAMPLE / CONTROL DATA
%{


rMx = randn(size(IMGsum)) .* .2;

rN = 30;
r1 = randn(rN) .* .5 + 4;
r2 = r1/2;
r3 = r1*2;

rMx(rN*3:rN*4-1,rN*3:rN*4-1) = r1;
rMx(rN*6:rN*7-1,rN*6:rN*7-1) = r2;
rMx(rN*9:rN*10-1,rN*9:rN*10-1) = r3;

    Mask = MaskMaker(1.2, 25, .4, .05);
rMx = convn(rMx,Mask,'same');



maxrMx = max(max(rMx));
minrMx = min(min(rMx));


lintrans = @(x,a,b,c,d) (c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)));

    for nn = 1:numel(rMx)

        x = rMx(nn);

        if maxrMx > 1
            rMx(nn) = lintrans(x,minrMx,maxrMx,0,1);
        else
            rMx(nn) = lintrans(x,minrMx,1,0,1);
        end


    end


figure(76)
imagesc(rMx)

IMGsum = rMx;

%}

%% -- NORMALIZE DATA TO RANGE: [0 <= DATA <= 1]

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Normalizing data to range [0 <= DATA <= 1]');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

IMGsumOrig = IMGsum;

IMGs = IMGsum;

% IMGs = IMGs .* 2;

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




%% ----- 1st FIGURES (ASIDE FROM OPTIONAL Z-STACK ANIMATION)  -----

    %fig3 = figure(3); set(fig3,'OuterPosition',[100 100 1100 500])
    %hax1 = axes('Position',[.05 .05 .4 .9]);
    %hax2 = axes('Position',[.52 .05 .4 .9]);

    % GUIfh = figure(NeuroFIT_GUI);
    % figure(NeuroFIT_GUI);

    %cla(GUIfh.Children(1).Children(2))
    %cla(GUIfh.Children(1).Children(1))
    %hax1 = GUIfh.Children(1).Children(1);
    %hax2 = GUIfh.Children(1).Children(2);

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Displaying original vs normalized image');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------



        %axes(hax1)
    axes(GUIfh.Children(1).Children(1));
imagesc(IMGsumOrig)

        %axes(hax2)
    axes(GUIfh.Children(1).Children(2));
imagesc(IMGs)

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

    MASKTBLR = [pos1(2) (pos1(2)+pos1(4)) pos1(1) (pos1(1)+pos1(3))];

    % Background
    mask{1} = zeros(size(iDUB));
    mask{1}(MASKTBLR(1):MASKTBLR(2), MASKTBLR(3):MASKTBLR(4)) = 1;
    mask1 = mask{1};


%% CHECK THAT MASK(S) ARE CORRECT

%     set(figure(77),'OuterPosition',[100 400 800 400],'Color',[1,1,1])
%     HaxP = [.04 .07 .42 .85; .53 .07 .42 .85];
% 
%         hax = axes('Position',HaxP(1,:));
%     %imagesc(iIND.*double(mask{1}));
%     imagesc(iDUB.*mask{1});

        axes(GUIfh.Children(1).Children(1));
    imagesc(iDUB.*mask{1});
        pause(2)


%% -- GET MEAN OF BACKGROUND PIXELS & SUBTRACT FROM IMAGE

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('Taking average of background pixels');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

    f1BACKGROUND = iDUB .* mask1;
    meanBG = mean(f1BACKGROUND(f1BACKGROUND > 0));

    meanALL = mean(iDUB(:));

    iDUB = iDUB - meanBG;
    %iDUB = iDUB - (meanBG / (meanBG+meanALL) * meanBG);
    iDUB(iDUB <= 0) = 0;



%% -- MESH SURFACE PLOT

    %GUIfh = figure(NeuroFIT_GUI);

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('displaying image intensity peaks as 3D mesh');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------


        axes(GUIfh.Children(1).Children(2));

    mesh(iDUB)
        view([186 24])

    pause(2)

%% check vague SNR for masks.

%     %-------
%     fig43 = figure(43);
%     set(fig43,'OuterPosition',[400  400  1200  600],'Color',[1,1,1])
%     %-------

     %axes(GUIfh.Children(1).Children(1));

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('plotting histogram of pixel intensities');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------


  hist((GUIfh.Children(1).Children(1)),iDUB(:),80);
    %xlim([.05 .9])
    %xlim([10 200])
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

    % select subset of pixels to use in terms of their brightness rank
	% this is up to the users discretion - and can be set in the dialogue box
	% the dialogue prompt has default vaules set to assess pixels that are 
	% within the range of X%-99.99% brightest


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
    spf0=sprintf('MEAN: % 4.4g ',mt);
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    spf0=sprintf('MEDIAN: % 4.4g ',mdn);
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    spf0=sprintf('STDEV: % 4.4g ',st);
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    pause(2)
    %----------------------------


%% -- Boxplot & Histogram

    %axes(GUIfh.Children(1).Children(2));

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('creating boxplot and histogram of analyzed pixels');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------


    boxplot((GUIfh.Children(1).Children(2)),Hi2LoIncPixArray ...
	,'notch','on' ...
	,'whisker',1 ...
	,'widths',.8 ...
	,'factorgap',[0] ...
	,'medianstyle','target');
	%set(gca,'XTickLabel',{' '},'Position',[.04 .05 .25 .9])
    pause(2)

    % axes(GUIfh.Children(1).Children(1));
hist((GUIfh.Children(1).Children(1)),Hi2LoIncPixArray(:),100);
		pause(2)

%% -- VIEW ORIGINAL, SNR MASK, AND TARGET IMAGES

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('showing original, snr mask, and target images');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

    %-------
    fig23 = figure(23);
    set(fig23,'OuterPosition',[10  50  1600  600],'Color',[1,1,1])
    %-------

        HaxOrig = axes('Position',   [.03 .05 .3 .9]);
    imagesc(IMGsumOrig);
        %axis(HaxOrig,'image')

        HaxThresh = axes('Position', [.36 .05 .3 .9]);
    imagesc(iDUB);
        %axis(HaxThresh,'image')

        HaxInc = axes('Position',    [.69 .05 .3 .9]);
    imagesc(IncRawMxP);
        %axis(HaxInc,'image')

    colormap('bone');

%% --- SET COLORMAP --

    cmaps = cmaplist;
%     [smap,butn] = listdlg('PromptString','Select a Colormap:',...
%                     'SelectionMode','single','ListSize',[200 200],...
%                     'ListString',cmaps);

   smap = 1;
    if smap == 1
        cmap = customcmap(1);
    else
        cmap = cmaps{smap};
    end
    colormap(cmap)
    pause(2)

    close gcf


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

            %fig25 = figure(25);
            %set(fig25,'OuterPosition',[400  100  1000  800],'Color',[1,1,1])

            axes(GUIfh.Children(1).Children(2))
            imagesc(fIMG)
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




cla(GUIfh.Children(1).Children(2))
cla(GUIfh.Children(1).Children(1))
hax1 = GUIfh.Children(1).Children(1);
hax2 = GUIfh.Children(1).Children(2);

%% --- PLOT AUTOTRACED EDGES ---


    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('plotting cell autotraces');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

        % axes(GUIfh.Children(1).Children(2))
        axes(hax1)
    imagesc(fIMG)
        colormap(hax1,hot)
        title('Convolution Processed Image')
        hold on
        

        for k=1:numel(fbIMG)
            plot(hax1,fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
        end
        hold on


pause(2)

        
        %axes(GUIfh.Children(1).Children(2))
        axes(hax2)
    imagesc(IMGsumOrig)
        colormap(hax2,hot)
        title('Original Image Sum')
        hold on

    for k=1:numel(fbIMG)
        plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w','LineWidth',2);
    end
        hold on

pause(2)

%%




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

    imgname = 'NFIT_IMG_512px_RGB_16bit.tif';

    [f2,map] = imread(imgname);				% import image
    infoto = imfinfo(imgname);              % image meta info
    iDUBr = im2double(f2);              % colormap is Mx-by-1 array of class double.
    iDUBr = imresize(iDUBr, [512 NaN]);
    iDUBr(iDUBr > 1) = 1;


%% -- CONVERT SECOND IMAGE TO USABLE MATRIX

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('processing complementary image');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

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

        axes(GUIfh.Children(1).Children(1))
    imagesc(iDUBr)
        title('Comparison Image')


%% -- SHOW FIGURE WITH IMAGE COMPARISONS

%     %-------
%     fig36 = figure(36);
%     set(fig36,'OuterPosition',[50  200  1500  800],'Color',[1,1,1])
%     %--
%     axL = [.03 .05 .45 .85];
%     axR = [.52 .05 .45 .85];
%     %-------
% 
%         gax = axes('Position',axL);
%     imagesc(iDUBg)
%         ccmapG = customcmap(2);
%         colormap(gax,ccmapG)
%         title('Green Channel')
%         hold on
%         
%     
%         for k=1:numel(fbIMG)
%     plot(fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
%         end
%         
% 
%         rax = axes('Position',axR);
%     imagesc(iDUBr)
%         ccmapR = customcmap(1);
%         colormap(rax,ccmapR)
%         title('Red Channel')
%         hold on
% 
%     for k=1:numel(fbIMG)
%         plot(fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
%     end
% pause(2)



%%
cla(GUIfh.Children(1).Children(2))
cla(GUIfh.Children(1).Children(1))
%%

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('displaying complimantary images');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------


        axes(hax1)
    imagesc(iDUBg)
        ccmapG = customcmap(2);
        colormap(hax1,ccmapG)
        title('Green Channel Processed Image')
        hold on
        

        for k=1:numel(fbIMG)
    plot(hax1,fbIMG{k}(:,2),fbIMG{k}(:,1),'g','LineWidth',2);
        end
        hold on


pause(.5)

        axes(hax2)
    imagesc(iDUBr)
        ccmapR = customcmap(1);
        colormap(hax2,ccmapR)
        title('Red Channel Comparison Image')
        hold on

    for k=1:numel(fbIMG)
        plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w','LineWidth',2);
    end
        hold on

pause(2)



%%
cla(GUIfh.Children(1).Children(2))
cla(GUIfh.Children(1).Children(1))
%%

    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('creating binary image');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

        % Create black and white (BW) binary image
    fIMG_BW = im2bw(fIMG_filled, .5);


    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('assigning labels to discrete objects in binary image');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

        % Returns labels for connected components of a binary image
    [L, NUM] = bwlabeln(fIMG_BW);  

        axes(hax1)
    imagesc(fIMG_BW);

        axes(hax2)
    imagesc(L);

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




% %%
% fontSize = 14;	% Used to control size of "blob number" labels put atop the image.
% labelShiftX = -7;	% Used to align the labels in the centers of the coins.
% blobECD = zeros(1, numberOfBlobs);
% % Print header line in the command window.
% fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
% % Loop over all blobs printing their measurements to the command window.
% for k = 1 : numberOfBlobs           % Loop through all blobs.
% 	% Find the mean of each blob.  (R2008a has a better way where you can pass the original image
% 	% directly into regionprops.  The way below works for all versions including earlier versions.)
%     thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
%     meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
% 	meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
% 	
% 	blobArea = blobMeasurements(k).Area;		% Get area.
% 	blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
% 	blobCentroid = blobMeasurements(k).Centroid;		% Get centroid.
% 	blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
%     fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));
% 	% Put the "blob number" labels on the "boundaries" grayscale image.
% 	text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold');
% end
% 
% %%


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

%%

cla(GUIfh.Children(1).Children(1))
cla(GUIfh.Children(1).Children(2))
% hax1 = GUIfh.Children(1).Children(1);
% hax2 = GUIfh.Children(1).Children(2);
%%


% %-------
% fig24 = figure(24);
% pos1 = [100  250  1200  800];
% set(fig24,'OuterPosition',pos1,'Color',[1,1,1])
	%HaxBP = axes('Position',[.04 .05 .25 .9]);
    %HaxHis = axes('Position',[.35 .05 .6 .9]);
% %-------


    %--- PRINT MESSAGE TO CON ---
    spf0=sprintf('plotting bar graph of pixel intensities for traced cells');
    [ft,spf2,spf3,spf4]=upcon(ft,spf0,spf2,spf3,spf4);
    %----------------------------

	
    axes(hax1), hold off
bar(PVs)

%     % Boxplot
%     axes(hax1), hold off
% hbp = boxplot(PVs ...
% 	,'notch','on' ...
% 	,'whisker',1 ...
% 	,'widths',.8 ...
% 	,'factorgap',[0] ...
% 	,'medianstyle','target');
% 	%set(HaxBP,'XTickLabel',{' '},'Position',[.04 .05 .25 .9],'YLIM',[0 1]);

	
    axes(hax2), hold off
imagesc(PLBW);


pause(1)

finalIMG = iDUB;
finalIMG(~PLBW) = 0;

        axes(hax2)
    imagesc(finalIMG)
        ccmapG = customcmap(2);
        colormap(hax2,ccmapG)
        title('Final Cell Trace')
        hold on
        

        for k=1:numel(fbIMG)
    plot(hax2,fbIMG{k}(:,2),fbIMG{k}(:,1),'w:','LineWidth',1.5);
        end
        hold on


pause(2)

%%


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
