function varargout = ImageMRIDat_DVD(varargin)
% IMAGEMRIDAT_DVD MATLAB code for ImageMRIDat_DVD.fig
%      IMAGEMRIDAT_DVD, by itself, creates a new IMAGEMRIDAT_DVD or raises the existing
%      singleton*.
%
%      H = IMAGEMRIDAT_DVD returns the handle to a new IMAGEMRIDAT_DVD or the handle to
%      the existing singleton*.
%
%      IMAGEMRIDAT_DVD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEMRIDAT_DVD.M with the given input arguments.
%
%      IMAGEMRIDAT_DVD('Property','Value',...) creates a new IMAGEMRIDAT_DVD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageMRIDat_DVD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageMRIDat_DVD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageMRIDat_DVD

% Last Modified by GUIDE v2.5 10-Nov-2023 19:10:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageMRIDat_DVD_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageMRIDat_DVD_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% addpath(genpath('C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\pvmatlab'));
addpath(genpath('pvmatlab'));


% End initialization code - DO NOT EDIT


% --- Executes just before ImageMRIDat_DVD is made visible.
function ImageMRIDat_DVD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageMRIDat_DVD (see VARARGIN)

% Choose default command line output for ImageMRIDat_DVD
handles.output = hObject;
axis(handles.FullImage);
imshow(zeros(2,2));
% Update handles structure
handles.GUIDataAll = [];
handles.GUIStartingPoint = handles;
guidata(hObject, handles);



set(gcf, 'units', 'normalized', 'position', [0 0 0.9 0.7])

% UIWAIT makes ImageMRIDat_DVD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImageMRIDat_DVD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in DataFolder.
function DataFolder_Callback(hObject, eventdata, handles)
% hObject    handle to DataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %% Define direcotry path for data, create temporary directory and re-arrange folders from ParaVision 360
    [fidfile, foldpath] = uigetfile('.\*.*');
%     addpath(genpath(foldpath));

    handles = guidata(hObject);
    handles.DataFolderPath.String = join([foldpath, fidfile],'');
    handles.GUIDataAll.CSIgenpath = foldpath;
    guidata(hObject, handles);
    
    if ~isfolder(join([foldpath, '\tmp_img'],''))
        mkdir(join([foldpath, '\tmp_img'],''))
    end

    filess = dir(fullfile(join([foldpath, '\pdata\1'],'')));
    for k = 3:size(filess)
        if ~isfile(join([foldpath, filess(k).name],''))
            copyfile(join([foldpath, '\pdata\1\', filess(k).name],''), join([foldpath, '\', filess(k).name],''))
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Process data to get right structure
    % Create a RawDataObject, importing the test data
    [rawObj, numSlices, numReps, acqSizes, fidFile2, fidFile3, numPhases, frameObj, ...
        kdataObj, dat, datsiz, kdataObj2, FIDdat, imageObj, imageObj2] = processMagDat(foldpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Phase Correction
    if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], ''))

        [imageObj3, sizz, epc, ephci, siss, phv, pap, phma, dsz, ph1, pivotppm, pivot, handles] = processPhaseDat(foldpath, imageObj2, handles);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        
        save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhaseParams.mat'], ''),'ephci', 'phma', 'dsz', 'ph1', 'pivotppm', 'pivot')
        handles = guidata(hObject);
        handles.GUIDataAll.ephci = ephci;
        handles.GUIDataAll.phma = phma;
        handles.GUIDataAll.dsz = dsz;
        handles.GUIDataAll.ph1 = ph1;
        handles.GUIDataAll.pivotppm = pivotppm;
        handles.GUIDataAll.pivot = pivot;
        handles.GUIDataAll.kdataObj2 = kdataObj2;
        guidata(hObject, handles);

    else
        imageObj3=kdataObj2.reco('quadrature');
        FD = load(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], '')).FDat3;
        imageObj3.data = FD;

        pp = load(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhaseParams.mat'], ''));
        handles = guidata(hObject);
        handles.GUIDataAll.ephci = pp.ephci;
        handles.GUIDataAll.phma = pp.phma;

        handles.GUIDataAll.dsz = pp.dsz;
        handles.GUIDataAll.ph1 = pp.ph1;
        handles.GUIDataAll.pivotppm = pp.pivotppm;
        handles.GUIDataAll.pivot = pp.pivot;
        guidata(hObject, handles);
    
    end


    handles = guidata(hObject);
    FDat = imageObj.data;
    FDat2 = imageObj2.data;
    FDat3 = imageObj3.data;
    handles.GUIDataAll.FDat = FDat;
    handles.GUIDataAll.FDat2 = FDat2;
    handles.GUIDataAll.FDat3 = FDat3;
    handles.GUIDataAll.Method = rawObj.Method;
    handles.GUIDataAll.numReps = numReps;
    handles.GUIDataAll.FIDDat = kdataObj2;
    
    guidata(hObject, handles);
 

    if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData.mat'], ''))
        save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData.mat'], ''),'FIDdat', 'FDat3', 'FDat2', 'FDat')
    end
    if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], ''))
        save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], ''),'FDat3')
    end


%% Find SNR value and Save Data
    [SNR, mva] = compSNR(handles, FDat, FDat3);
    
    handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
    guidata(hObject, handles);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if ~isfile(join([foldpath, 'tmp_img\CSIProcessedData.mat'], ''))
        Magnitude = FDat;
        Complex = imageObj2.data;
        AutoPhased = imageObj3.data;
        save(join([foldpath, 'tmp_img\CSIProcessedData.mat'], ''), "Magnitude", "Complex", 'AutoPhased');
    end
    

    % This save each spectra in a different CSV file in case someone wants
    % to work with the data without having to deal with matlab, but it
    % takes a really long time for big CSIs. Perhaps put everything into
    % one same CSV with a header indicating the row and column??????

%     if ~isfolder(join([foldpath, 'tmp_img\ProcessedDataCSV'], ''))
%         mkdir(join([foldpath, 'tmp_img\ProcessedDataCSV'], ''))
%         datsiz2 = size(FDat);
%         for i = 1:datsiz2(2)
%             for j = 1:datsiz2(3)
%                 for k = 1:numReps
%                     T = array2table([real(imageObj2.data(:,i,j,1,1,1,k)), imag(imageObj2.data(:,i,j,1,1,1,k)), FDat(:,i,j,1,1,1,k)]);
%                     T.Properties.VariableNames(1:3) = {'Real','Imaginary','Magnitude'};
%                     writetable(T,join([foldpath, 'tmp_img\ProcessedDataCSV\CSIData_X',num2str(i),'_Y',num2str(j),'_T',num2str(k),'.csv'], ''))
%                 end
%             end
%         end
%     end

    



%% Generation of Grid

axes(handles.MainPlot)  
hold on

if rawObj.Method.PVM_Matrix(1) == rawObj.Method.PVM_Matrix(2) % Squared grids
    xAxisPix=0:rawObj.Method.PVM_Matrix(1);
    yAxisPix=0:rawObj.Method.PVM_Matrix(2);
    plot([xAxisPix;xAxisPix],repmat([yAxisPix(1);yAxisPix(end)],1,length(xAxisPix)),'Color', '#EDB120','LineWidth',1.5)
    plot(repmat([xAxisPix(1);xAxisPix(end)],1,length(xAxisPix)),[yAxisPix;yAxisPix],'Color', '#EDB120','LineWidth',1.5)
    handles.MainPlot.XLim = [xAxisPix(1) xAxisPix(end)];
    handles.MainPlot.YLim = [yAxisPix(1) yAxisPix(end)];
else % rectangular grids
    xAxisPix=0:rawObj.Method.PVM_Matrix(2);
    yAxisPix=0:rawObj.Method.PVM_Matrix(1);

    plot([xAxisPix;xAxisPix],repmat([yAxisPix(1);yAxisPix(end)],1,length(xAxisPix)),'Color', '#EDB120','LineWidth',1.5)
    plot(repmat([xAxisPix(1);xAxisPix(end)],1,length(yAxisPix)),[yAxisPix;yAxisPix],'Color', '#EDB120','LineWidth',1.5)
    handles.MainPlot.XLim = [xAxisPix(1) xAxisPix(end)];
    handles.MainPlot.YLim = [yAxisPix(1) yAxisPix(end)];

end


TPs = numReps;

handles = guidata(hObject);
handles.TimePoints.String = [1:TPs];
guidata(hObject, handles);

    
% Generate the actual CSI plot
genSpectraAll(FDat, handles.TimePoints.Value, foldpath)
answerRC = genSpectraAll_Phased(real(FDat3), handles.TimePoints.Value, foldpath);

handles = guidata(hObject);
handles.GUIDataAll.answerRC = answerRC;
guidata(hObject, handles);


if handles.togglebutton7.Value == 0
    im = imread(join([foldpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
else
    im = imread(join([foldpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'_Phased.png'], ''));
end
axes(handles.SpectraPan)

imm=imshow(im);
imm.Parent.XLim = [1-0.5   800+0.5];
imm.Parent.YLim = [1-0.5   800+0.5];
imm.XData = [1   800];
imm.YData = [1   800];



    try   
        % get the visu parameters for the grid
        [FOVSizeCSI, FOVPosCSI] = getVisuPars(handles, foldpath);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        if handles.FullImage.XLim(2)>5
            % Extract Field of View Coordenates for CSI and Proton Images
            % THIS IS UNDER TESTING -- DAVID
            [posIM, posCS, posCSor, posIMor, xb, yb] = getFOVCors(handles, FOVPosCSI, FOVSizeCSI);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            handles = guidata(hObject);
            handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
            handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
            guidata(hObject, handles);
            % Check if it is using the turbo localiser so we can flip the x
            % axis
    
            if (posCSor(1)<0 && posIMor(1)>0) || (posCSor(1)>0 && posIMor(1)<0) %turbo == 1%turbo == 1
%             xb = impix(1) - xb;

                Y = flip(handles.GUIDataAll.MRIImage,2);
                handles.GUIDataAll.MRIImage = Y;
                guidata(hObject, handles);
                axes(handles.FullImage)
                imshow(Y,[])
    
            end
    
            if (posCSor(2)<0 && posIMor(2)>0) || (posCSor(2)>0 && posIMor(2)<0)
                Y = flip(handles.GUIDataAll.MRIImage,1);
                handles.GUIDataAll.MRIImage = Y;
                guidata(hObject, handles);
                axes(handles.FullImage)
                imshow(Y,[])
            end

            axes(handles.FullImage)
            hold on 
            patch(xb,yb, 'g', 'EdgeColor', 'none','FaceAlpha',.4);

        end

        handles = guidata(hObject);
        handles.GUIDataAll.CSIPos = FOVPosCSI;
        handles.GUIDataAll.CSISize = FOVSizeCSI;
        handles.GUIDataAll.CSIVox = size(FDat);
        guidata(hObject, handles);
        
        
        if isfield(handles.GUIDataAll, 'MRIImage')
            if ~isempty(handles.GUIDataAll.MRIImage)
            MRIIm = handles.GUIDataAll.MRIImage;
            cutMRIIm = MRIIm(max([min(yb) 1]):max(yb), max([min(xb) 1]):max(xb));
            handles = guidata(hObject);
            handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
            handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
            guidata(hObject, handles);
            imwrite(mat2gray(cutMRIIm), join([foldpath, '\tmp_img\MRICut.png'], ''));    

            mim = size(FDat);

            axes(handles.MainPlot)
            fd = imshow(mat2gray(cutMRIIm));
            fd.XData = [0 mim(2)];
            fd.YData = [0 mim(3)];
            fd.Parent.XLim = [0,mim(2)];
            fd.Parent.YLim = [0,mim(3)];
            h = get(handles.MainPlot,'Children');
            set(handles.MainPlot,'Children',[h(2:end); h(1)])

            end
        end
        
    catch
        f = warndlg('No visu_pars file in the same folder than the fid file! Please, fix this if you want to overlay CSI and MRI Image!','Warning');
    end

% Compute the ppm axix for the plots
try
    bw = imageObj2.Method.PVM_SpecSW(1);
catch
    bw = imageObj2.Method.SpecBandPpm(1);
end
bwc = imageObj2.Method.PVM_FrqWorkPpm(1);

   
ppms = flip(linspace(bwc-bw/2, bwc+bw/2, datsiz(1)));

if ~isfile(join([foldpath, 'tmp_img\ppms.mat'], ''))
    save(join([foldpath, 'tmp_img\ppms.mat'], ''), "ppms");
end
handles.GUIDataAll.ppms = ppms;
guidata(hObject, handles);


%% Selection of voxels and ploting of spectra
try
while true
    
    pause(0.5)
    stop_state = get(handles.ClearGUI, 'Value');
    if stop_state
        handles = guidata(hObject);
        handles.ClearGUI.Value = 0;
        guidata(hObject, handles);
        break;
    end

    [x,y, but] = ginput(1);
    
    x2=floor(x)+1;
    y2=floor(y)+1;


    if but == 3
        handles.PPM.String = x;
            [minValue,closestIndex] = min(abs(ppms-x));
            handles.Intens.String = handles.GUIDataAll.currDat(closestIndex);
    end

    xa = [x2-1 x2 x2 x2-1];
    
    try
        fs.FaceColor='none';
    catch
    end

    h = get(handles.MainPlot,'Children');
    if isempty(findobj(h, 'type', 'Image'))
        y2 = (handles.GUIDataAll.CSIVox(3) - y2)+1;
        ya = [handles.GUIDataAll.CSIVox(3)-y2 handles.GUIDataAll.CSIVox(3)-y2 handles.GUIDataAll.CSIVox(3)-y2+1 handles.GUIDataAll.CSIVox(3)-y2+1];
    else
        ya = [y2-1 y2-1 y2 y2];
    end

    
    axes(handles.MainPlot)
    fs=patch(xa,ya, 'g', 'EdgeColor', 'none','FaceAlpha',.2);
    
    if (x2 <= 0) || (x2 >= handles.GUIDataAll.CSIVox(2)+1) || (y2<=0) || (y2 >= handles.GUIDataAll.CSIVox(3)+1)
        handles = guidata(hObject);
        try 
            x2 = handles.GUIDataAll.x2;
            y2 = handles.GUIDataAll.y2;
            handles.XYPos.String = join(['X: ',string(x2),', Y: ',string(y2),''], '');
        catch
            handles.XYPos.String = join(['X: ?, Y: ?'], '');
        end
        guidata(hObject, handles);
    else
        
        handles = guidata(hObject);
        handles.XYPos.String = join(['X: ',string(x2),', Y: ',string(y2),''], '');

        handles.GUIDataAll.x2 = x2;
        handles.GUIDataAll.y2 = y2;
        guidata(hObject, handles);
        
    end

    guidata(hObject, handles);
    rd = real(handles.GUIDataAll.FDat2(:,:,:,1,1,1,handles.TimePoints.Value));
    id = imag(handles.GUIDataAll.FDat2(:,:,:,1,1,1,handles.TimePoints.Value));
    if handles.togglebutton7.Value == 0
        fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,handles.TimePoints.Value);
    else
        fd = real(handles.GUIDataAll.FDat3(:,:,:,1,1,1,handles.TimePoints.Value));
    end
    

    dimsdat = size(imageObj.data);

    % Actual Ploting

    try 
        if ~isempty(handles.axes2.Children)
            set(handles.axes2.Children, 'YData', real(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))
            set(handles.axes2.Children, 'XData', ppms)
            set(handles.axes2, 'XLim', [min(ppms) max(ppms)])
            set(handles.axes2, 'YLim', [min(rd(:)) max(rd(:))])
            set(handles.axes2, 'XDir', 'reverse')
        else
            axes(handles.axes2)  
            cla reset
            plot(ppms, real(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
            ylim([min(rd(:)) max(rd(:))])
            xlim([min(ppms) max(ppms)])
            ax = gca;
            ax.XDir = 'reverse';
    %         xlabel('ppm')
            legend('Real')
        end
        
        if ~isempty(handles.axes3.Children)
            set(handles.axes3.Children, 'YData', imag(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))
            set(handles.axes3.Children, 'XData', ppms)
            set(handles.axes3, 'XLim', [min(ppms) max(ppms)])
            set(handles.axes3, 'YLim', [min(id(:)) max(id(:))])
            set(handles.axes3.XLabel, 'String', 'ppm')
            set(handles.axes3, 'XDir', 'reverse')
        else
            axes(handles.axes3)
            cla reset
            plot(ppms,imag(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
            ylim([min(id(:)) max(id(:))])
            xlim([min(ppms) max(ppms)])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')
            legend('Imag')
        end
        
        if isnan(str2double(handles.SFac.String))

        if ~isempty(handles.axes4.Children)

            if handles.togglebutton7.Value == 0 
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))
                set(handles.axes4.Children, 'XData', ppms+str2double(handles.CSC.String))
                set(handles.axes4.Legend, "String", 'Magn')
            else
                set(handles.axes4.Children, 'YData', real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))
                set(handles.axes4.Children, 'XData', ppms+str2double(handles.CSC.String))
                set(handles.axes4.Legend, 'String', 'Phs. Corr.')
            end

            set(handles.axes4, 'XLim', [min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])

        else

            axes(handles.axes4)
            cla reset
            if handles.togglebutton7.Value == 0 
                plot(ppms+str2double(handles.CSC.String),handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))), 'g')
                legend('Magn')
            else
                plot(ppms+str2double(handles.CSC.String),real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'g')
                legend('Phs. Corr.')
            end
            ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
            xlim([min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')
        end
        
            handles = guidata(hObject);
        
            if handles.togglebutton7.Value == 0 
                handles.GUIDataAll.currDat = handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            else
                handles.GUIDataAll.currDat = handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            end
            guidata(hObject, handles);

        
        else

        if ~isempty(handles.axes4.Children)

            if handles.togglebutton7.Value == 0 
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms+str2double(handles.CSC.String))
                set(handles.axes4.Legend, "String", 'Magn')
            else
                set(handles.axes4.Children, 'YData', real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms+str2double(handles.CSC.String))
                set(handles.axes4.Legend, 'String', 'Phs. Corr.')
            end

            set(handles.axes4, 'XLim', [min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])

        else

            axes(handles.axes4)
            cla reset
            if handles.togglebutton7.Value == 0  
                plot(ppms+str2double(handles.CSC.String),handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
                legend('Magn')
            else
                plot(ppms+str2double(handles.CSC.String),real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String), 'g')
                legend('Phs. Corr.')
            end
            ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
            xlim([min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')

        end
            
            handles = guidata(hObject);
       
            if handles.togglebutton7.Value == 0 
                handles.GUIDataAll.currDat = handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            else
                handles.GUIDataAll.currDat = handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            end
             guidata(hObject, handles);
        end
        
    
    catch
        
    end

    
end

catch
    disp('GUI closed')
end




function DataFolderPath_Callback(hObject, eventdata, handles)
% hObject    handle to DataFolderPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataFolderPath as text
%        str2double(get(hObject,'String')) returns contents of DataFolderPath as a double


% --- Executes during object creation, after setting all properties.
function DataFolderPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataFolderPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TimePoints.
function TimePoints_Callback(hObject, eventdata, handles)
% hObject    handle to TimePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    genSpectraAll(handles.GUIDataAll.FDat, handles.TimePoints.Value, handles.GUIDataAll.CSIgenpath)
    answerRC = genSpectraAll_Phased(real(handles.GUIDataAll.FDat3), handles.TimePoints.Value, handles.GUIDataAll.CSIgenpath);

    if handles.togglebutton7.Value == 0
        im = imread(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
    else
        im = imread(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'_Phased.png'], ''));
    end
    axes(handles.SpectraPan)
    imshow(im)

    % Find SNR value
    [SNR, mva] = compSNR(handles, handles.GUIDataAll.FDat, handles.GUIDataAll.FDat3);
    
    handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
    guidata(hObject, handles);

    handles = guidata(hObject);

    handles.GUIDataAll.answerRC = answerRC;

    handles.togglebutton1.Value = 0;
    togglebutton1_Callback(hObject, [], handles)

    handles.togglebutton4.Value = 0;
    togglebutton4_Callback(hObject, [], handles)

    guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns TimePoints contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TimePoints


% --- Executes during object creation, after setting all properties.
function TimePoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushTest.
function pushTest_Callback(hObject, eventdata, handles)
% hObject    handle to pushTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    [fidfile, foldpath] = uigetfile('.\*.*');
%     addpath(genpath(foldpath));
    
    
    if ~strcmp(lower(fidfile(end-3:end)),'.dcm')
        

        filess = dir(fullfile(join([foldpath, '\pdata\1'],'')));
        for k = 3:size(filess)
            if ~isfile(join([foldpath, filess(k).name],''))
                copyfile(join([foldpath, '\pdata\1\', filess(k).name],''), join([foldpath, '\', filess(k).name],''))
            end
        end

        kdataObj = CKDataObjectIMAGE(foldpath);


        if strcmp(kdataObj.Acqp.ACQ_protocol_name, 'EPI')
            f = warndlg('Sorry, but processing of EPI images is still not supported. We are working on it');
            waitfor(f)

            % Need to figure out why the processing gives folding
            % artifacts, while when done in ParaVision it doesnt. 
            
            Y = zeros(size(kdataObj.data));
            handles = guidata(hObject);
            handles.edit2.String = join([foldpath, fidfile],'');
            axes(handles.FullImage)
            imshow(Y,[])
           
            handles.GUIDataAll.MRIImage = Y;
            guidata(hObject, handles);

        else
            f = warndlg('DICOM file NOT selected! Image will be reconstructed from FID file. This is NOT the preferred option! If more than one plane is acquired, only the last one acquired will be displayed!');
            waitfor(f)
    
            kdataObj = kdataObj.readReco;
            imageObj=kdataObj.reco('all', 'image');
    
    
            handles = guidata(hObject);
            handles.edit2.String = join([foldpath, fidfile],'');
            Y = mat2gray(imageObj.data(:,:,1,1,1,end)');
            axes(handles.FullImage)
            imshow(Y,[])
           
            handles.GUIDataAll.MRIImage = Y;
            guidata(hObject, handles);
        end

    else

        handles = guidata(hObject);
        handles.edit2.String = join([foldpath, fidfile],'');
        Y = dicomread(join([foldpath, fidfile],''));
        axes(handles.FullImage)
        imshow(Y,[])
        
        handles.GUIDataAll.MRIImage = Y;
        guidata(hObject, handles);
    end

    
    try % Try to oppen the visu_pars file for the localiser to get inforamtion of the FOV size
        
        [FOVSizeIm, FOVPosIm, turbo] = getVisuParsH(foldpath);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    catch % Assume a FOV size and display a warning window telling it to the user
        
        FOVSizeIm = [60 60];
        FOVPosIm = [0 0 0];
        turbo = 0;

        f = warndlg('Could not find the visu_pars file, so there is no information about the FOV size. The code will assume 60 x 60 mm','Warning');
        
    end
    handles = guidata(hObject);
    handles.GUIDataAll.MRIImage = Y;
    handles.GUIDataAll.MRIImageSize = FOVSizeIm;
    handles.GUIDataAll.MRIImagePos = FOVPosIm(end-2:end);
    handles.GUIDataAll.turbo = turbo;
    guidata(hObject, handles);

    


    if isfield(handles.GUIDataAll, 'CSIPos')
        if ~isempty(handles.GUIDataAll.CSIPos)

            [impix, FOVSizeCSI, FOVPosCSI, posIM, posCS, posCSor, posIMor] = getFOVCorsH(handles);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (posCSor(1)<0 && posIMor(1)>0) || (posCSor(1)>0 && posIMor(1)<0) %turbo == 1

            Y = flip(dicomread(join([foldpath, fidfile],'')),2);
            handles.GUIDataAll.MRIImage = Y;
            guidata(hObject, handles);

        end

        if (posCSor(2)<0 && posIMor(2)>0) || (posCSor(2)>0 && posIMor(2)<0)
            Y = flip(dicomread(join([foldpath, fidfile],'')),1);
            handles.GUIDataAll.MRIImage = Y;
            guidata(hObject, handles);
        end

        axes(handles.FullImage)
        imshow(Y,[])

        PosGen = posIM-posCS;
        
        handles.GUIDataAll.PosGen = PosGen;

        impix = size(handles.GUIDataAll.MRIImage);

        if impix(1) == impix(2)
            xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
            yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
        else
            yb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
            xb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
        end

        % Check if it is using the turbo localiser so we can flip the x
        % axis

        axes(handles.FullImage)
        hold on 
        patch(xb,yb,'g', 'EdgeColor', 'none','FaceAlpha',.4);
        
        
        MRIIm = handles.GUIDataAll.MRIImage;
        cutMRIIm = MRIIm(max([min(yb) 1]):max(yb), max([min(xb) 1]):max(xb));
        imwrite(mat2gray(cutMRIIm), join([handles.GUIDataAll.CSIgenpath, '\tmp_img\MRICut.png'], '')); 
        handles = guidata(hObject);
        handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
        handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
        guidata(hObject, handles);
        mim = handles.GUIDataAll.CSIVox;

        axes(handles.MainPlot)
        fd = imshow(mat2gray(cutMRIIm));
        fd.XData = [0 mim(2)];
        fd.YData = [0 mim(3)];
        fd.Parent.XLim = [0,mim(2)];
        fd.Parent.YLim = [0,mim(3)];
        h = get(handles.MainPlot,'Children');
        set(handles.MainPlot,'Children',[h(2:end); h(1)])
        end
    end
    
    handles = guidata(hObject);
    handles.togglebutton1.Value = 0;
    togglebutton1_Callback(hObject, [], handles)

    handles.GUIDataAll.MRIImage = Y;
    
    guidata(hObject, handles);


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

if handles.togglebutton1.Value == 1
    handles = guidata(hObject);
    handles.togglebutton1.BackgroundColor = [0    1    0];
    guidata(hObject, handles);        

    if ~isempty(handles.DataFolderPath.String) && isempty(handles.edit2.String)

        foldpath = handles.GUIDataAll.CSIgenpath;

        if handles.togglebutton7.Value == 0
            imOver = imread(join([foldpath, '/tmp_img/SpectraTimePoint_',num2str(handles.TimePoints.Value),'.png'], ''));
        else

            imOver = imread(join([foldpath, '/tmp_img/SpectraTimePoint_',num2str(handles.TimePoints.Value),'_Phased.png'], ''));
        end

        mim = handles.GUIDataAll.CSIVox;
 
        axes(handles.MainPlot)
        h = get(handles.MainPlot,'Children');
        delete(findobj(h, 'type', 'Patch'));
        delete(findobj(h, 'type', 'Image'));

        handles = guidata(hObject);
        fd = imshow(imOver);
        fd.XData = [0 mim(2)];
        fd.YData = [0 mim(3)];
        fd.Parent.XLim = [0,mim(2)];
        fd.Parent.YLim = [0,mim(3)];
        h = get(handles.MainPlot,'Children');
        set(handles.MainPlot,'Children',[h(2:end); h(1)])
        guidata(hObject, handles);
% 
% 
%     else
    elseif ~isempty(handles.DataFolderPath.String) && ~isempty(handles.edit2.String)
        genOverlayPythonScript(handles.GUIDataAll.CSIgenpath, handles.TimePoints.Value);
        genOverlayPythonScript_Phased(handles.GUIDataAll.CSIgenpath, handles.TimePoints.Value);

        foldpath = handles.GUIDataAll.CSIgenpath;
        try
            answerRC = handles.GUIDataAll.answerRC;
        catch
            answerRC = [];
        end
        if ~isfile(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'],''))            
            system(['python "',foldpath,'/tmp_img/',num2str(handles.TimePoints.Value),'-MRICSIOverlayGen.py"']);
            system(['python "',foldpath,'/tmp_img/',num2str(handles.TimePoints.Value),'-MRICSIOverlayGen_Phased.py"']);
        end
        if isfile(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'],'')) && strcmp(answerRC, 'Yes')
            system(['python "',foldpath,'/tmp_img/',num2str(handles.TimePoints.Value),'-MRICSIOverlayGen_Phased.py"']);
            handles = guidata(hObject);
            handles.GUIDataAll.answerRC = [];
            guidata(hObject, handles);
        end
        
        impix = size(handles.GUIDataAll.MRIImage);
        FOVSizeCSI = handles.GUIDataAll.CSISize;
        FOVPosCSI = handles.GUIDataAll.CSIPos;
        FOVPosIm = handles.GUIDataAll.MRIImagePos;
        foldpath = handles.GUIDataAll.CSIgenpath;

        if handles.togglebutton7.Value == 0
            imOver = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'], ''));
        else
            imOver = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'_Phased.png'], ''));
        end
        

        if ~isfile(join([foldpath, '\tmp_img\MRICSIOverlay_FullImage.png'], ''))
        % Generate overlay with full image
            if handles.togglebutton7.Value == 0
                imOver2 = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'], ''));
            else
                imOver2 = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'_Phased.png'], ''));
            end
            imwrite(mat2gray(handles.GUIDataAll.MRIImage), join([foldpath, '\tmp_img\MRIFullImage.png'], '')); 
    
            % To make sure you get the grid edges in yellow instead of just
            % black
            imOver2(1,:,1) = 237; imOver2(1,:,2) = 177; imOver2(1,:,3) = 32;
            imOver2(:,1,1) = 237; imOver2(:,1,2) = 177; imOver2(:,1,3) = 32;
            imOver2(:,end,1) = 237; imOver2(:,end,2) = 177; imOver2(:,end,3) = 32;
    
            imcut = handles.GUIDataAll.MRIImage(round(handles.GUIDataAll.xIndCut), round(handles.GUIDataAll.yIndCut));
            imorg = imresize(cat(3, im2uint16(mat2gray(handles.GUIDataAll.MRIImage)), im2uint16(mat2gray(handles.GUIDataAll.MRIImage)), ...
                    im2uint16(mat2gray(handles.GUIDataAll.MRIImage))), size(handles.GUIDataAll.MRIImage)*8);
    
            imorg(round(handles.GUIDataAll.xIndCut(1))*8-7:round(handles.GUIDataAll.xIndCut(end))*8, ...
                    round(handles.GUIDataAll.yIndCut(1))*8-7:round(handles.GUIDataAll.yIndCut(end))*8,:) = imresize(im2uint16(imOver2), size(imcut)*8);
            
            imwrite(imorg, join([foldpath, '\tmp_img\MRICSIOverlay_FullImage.png'], ''));
        end


%         PosGen = abs(FOVPosCSI - FOVPosIm);

%         xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);
%         yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.GUIDataAll.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.GUIDataAll.MRIImageSize(2);

        mim = handles.GUIDataAll.CSIVox;

        axes(handles.MainPlot)
        h = get(handles.MainPlot,'Children');
        delete(findobj(h, 'type', 'Patch'));
        delete(findobj(h, 'type', 'Image'));

        handles = guidata(hObject);
        fd = imshow(imOver);
        fd.XData = [0 mim(2)];
        fd.YData = [0 mim(3)];
        fd.Parent.XLim = [0,mim(2)];
        fd.Parent.YLim = [0,mim(3)];
        h = get(handles.MainPlot,'Children');
        set(handles.MainPlot,'Children',[h(2:end); h(1)])
        guidata(hObject, handles);

    else
%         f = warndlg('To overlay the CSI image to the MRI DICOM image you need to load both files! :S','Warning');
    end
else
    handles = guidata(hObject);
    handles.togglebutton1.BackgroundColor = [0.9400    0.9400    0.9400];
    guidata(hObject, handles);

    if ~isempty(handles.DataFolderPath.String) && isempty(handles.edit2.String)

        handles = guidata(hObject);
        axes(handles.MainPlot)
        h = get(handles.MainPlot,'Children');
        delete(findobj(h, 'type', 'Patch'));
        delete(findobj(h, 'type', 'Image'));
        guidata(hObject, handles);
%         plot([0, 16],[0, 16])

% 
%     else
    elseif ~isempty(handles.DataFolderPath.String) && ~isempty(handles.edit2.String)
        
        [impix, FOVSizeCSI, FOVPosCSI, FOVPosIm, foldpath, turbo, posIM, posCS, PosGen, xb, yb] = getFOV(handles);


        MRIIm = handles.GUIDataAll.MRIImage;
        cutMRIIm = MRIIm(max([min(yb) 1]):max(yb), max([min(xb) 1]):max(xb));
        imwrite(mat2gray(cutMRIIm), join([foldpath, '\tmp_img\MRICut.png'], '')); 

        handles = guidata(hObject);
        handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
        handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
        guidata(hObject, handles);
        mim = handles.GUIDataAll.CSIVox;

        handles = guidata(hObject);
        axes(handles.MainPlot)
        h = get(handles.MainPlot,'Children');
        delete(findobj(h, 'type', 'Patch'));
        delete(findobj(h, 'type', 'Image'));

        fd = imshow(mat2gray(cutMRIIm));
        fd.XData = [0 mim(2)];
        fd.YData = [0 mim(3)];
        fd.Parent.XLim = [0,mim(2)];
        fd.Parent.YLim = [0,mim(3)];
        h = get(handles.MainPlot,'Children');
        set(handles.MainPlot,'Children',[h(2:end); h(1)])
        guidata(hObject, handles);


    else
%         f = warndlg('To overlay the CSI image to the MRI DICOM image you need to load both files! :S','Warning');
    end

end


% --- Executes on button press in ClearGUI.
function ClearGUI_Callback(hObject, eventdata, handles)
% hObject    handle to ClearGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guidata(hObject);
% 
% guidata(hObject, handles);


close(gcbf)
ImageMRIDat_DVD;

% handles = guidata(hObject);
% delete(handles.MainPlot.Children);
% handles.DataFolderPath.String = '';
% handles.edit2.String = '';
% handles.togglebutton1.Value = 0;
% handles.togglebutton1.BackgroundColor = [0.9400    0.9400    0.9400];
% delete(handles.FullImage.Children)
% delete(handles.SpectraPan.Children)
% delete(handles.axes2.Children)
% delete(handles.axes3.Children)
% delete(handles.axes4.Children)
% handles.GUIDataAll.CSIgenpath = [];
% handles.GUIDataAll.FDat = [];
% handles.GUIDataAll.CSIPos = [];
% handles.GUIDataAll.CSISize = [];
% handles.GUIDataAll.CSIVox = [];
% handles.GUIDataAll.MRIImage = [];
% handles.GUIDataAll.MRIImageSize = [];
% handles.GUIDataAll.MRIImagePos = [];
% handles.FullImage.XLim = [0 1];
% % handles.ClearGUI.Value = 0;
% 
% % handles = handles.GUIStartingPoint;
% 
% guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over ClearGUI.
function ClearGUI_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ClearGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.ClearGUI.Value = 0;
% delete(handles.MainPlot.Children);
% handles.DataFolderPath.String = '';
% handles.edit2.String = '';
% handles.togglebutton1.Value = 0;
% handles.togglebutton1.BackgroundColor = [0.9400    0.9400    0.9400];
% delete(handles.FullImage.Children)
% delete(handles.SpectraPan.Children)
% delete(handles.axes2.Children)
% delete(handles.axes3.Children)
% delete(handles.axes4.Children)
% guidata(hObject, handles);



function SFac_Callback(hObject, eventdata, handles)
% hObject    handle to SFac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SFac as text
%        str2double(get(hObject,'String')) returns contents of SFac as a double

if isnan(str2double(handles.SFac.String))
    f = warndlg('Please, enter a number or this is not going to work! :S','Warning');
else
    if isfield(handles.GUIDataAll, 'FDat')
        dimsdat = size(handles.GUIDataAll.FDat);
        if handles.togglebutton7.Value == 0
            fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        else
            fd = handles.GUIDataAll.FDat3(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        end
        x2 = handles.GUIDataAll.x2;
        y2 = handles.GUIDataAll.y2;
        ppms = handles.GUIDataAll.ppms;

        if ~isempty(handles.axes4.Children)

            if handles.togglebutton7.Value == 0 
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, "String", 'Magn')
            else
                set(handles.axes4.Children, 'YData', real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, 'String', 'Phs. Corr.')
            end

            set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
            set(handles.axes4, 'XLim', [min(ppms) max(ppms)])

        else

            axes(handles.axes4)
            cla reset
    %         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
            if handles.togglebutton7.Value == 0
                plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
                legend('Magn')
            else
                plot(ppms, real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String), 'g')
                legend('Phs. Corr.')
            end
            ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
            xlim([min(ppms) max(ppms)])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')
        end

    end
end



% --- Executes during object creation, after setting all properties.
function SFac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SFac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DelTmpDir.
function DelTmpDir_Callback(hObject, eventdata, handles)
% hObject    handle to DelTmpDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DelTmpDir

if isfield(handles.GUIDataAll, 'CSIgenpath')
    if isfolder(join([handles.GUIDataAll.CSIgenpath, '\tmp_img'],''))
        answer = questdlg('Are you sure you want to delete all temporary files? This will also close the GUI!', ...
	        'Decisions, decisions, decisions', ...
	        'Yes','No','No');
        if isempty(answer)
            answer = 'No';
        end
        if strcmp(answer, 'Yes')    
            delete(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\*'],''))
            closereq(); 
        end
    end
end

handles.DelTmpDir.Value = 0;


% --- Executes on button press in togglebutton4.
function togglebutton4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton4

if handles.togglebutton6.Value == 1
    handles = guidata(hObject);
    handles.togglebutton6.Value = 0;
    handles.togglebutton6.BackgroundColor = [0.9400    0.9400    0.9400];
    togglebutton6_Callback(hObject, [], handles);
    guidata(hObject, handles);
end


if handles.togglebutton4.Value == 1
    handles = guidata(hObject);
    handles.togglebutton4.BackgroundColor = [0    1    0];
    guidata(hObject, handles); 

    if isfield(handles, "GUIDataAll") && isfield(handles.GUIDataAll, "FDat")
%         genSpectraAllColorGrid(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String), handles.GUIDataAll.CSIgenpath);
        axes(handles.SpectraPan)
        genSpectraAllColorGrid_Colormap(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String(handles.TimePoints.Value,:)), handles.GUIDataAll.CSIgenpath, str2num(handles.edit8.String));

        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_ColorGrid_',string(handles.TimePoints.Value),'_Sat',string(handles.edit8.String),'.png'], ''));
        axes(handles.SpectraPan)

        rati = handles.GUIDataAll.CSISize(2:3)/max(handles.GUIDataAll.CSISize(2:3));

        imm=imshow(im);
        imm.Parent.XLim = [1-0.5   (800+0.5)*rati(1)];
        imm.Parent.YLim = [1-0.5   (800+0.5)*rati(2)];
        imm.XData = [1   800*rati(1)];
        imm.YData = [1   800*rati(2)];


    else
        f = warndlg('You need to load a CSI experiment first if you wanna generate the colored grid! :S','Warning');
        handles = guidata(hObject);
        handles.togglebutton4.Value = 0;
        handles.togglebutton4.BackgroundColor = [0.9400    0.9400    0.9400];
        guidata(hObject, handles);
    end

else
    handles = guidata(hObject);
    handles.togglebutton4.BackgroundColor = [0.9400    0.9400    0.9400];
    guidata(hObject, handles);

    if handles.togglebutton7.Value == 0
        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
    else
        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'_Phased.png'], ''));
    end
    axes(handles.SpectraPan)
    imm=imshow(im);
    imm.Parent.XLim = [1-0.5   800+0.5];
    imm.Parent.YLim = [1-0.5   800+0.5];
    imm.XData = [1   800];
    imm.YData = [1   800];

end



function CSC_Callback(hObject, eventdata, handles)
% hObject    handle to CSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CSC as text
%        str2double(get(hObject,'String')) returns contents of CSC as a double
if isnan(str2double(handles.CSC.String))
    f = warndlg('Please, enter a number or this is not going to work! :S','Warning');
else
    if isfield(handles.GUIDataAll, 'FDat')
        dimsdat = size(handles.GUIDataAll.FDat);
        if handles.togglebutton7.Value == 0
            fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        else
            fd = handles.GUIDataAll.FDat3(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        end
        x2 = handles.GUIDataAll.x2;
        y2 = handles.GUIDataAll.y2;
        ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);

        if ~isempty(handles.axes4.Children)

            if handles.togglebutton7.Value == 0 
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, "String", 'Magn')
            else
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, 'String', 'Phs. Corr.')
            end

            set(handles.axes4, 'XLim', [min(ppms) max(ppms)])
            set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])

        else

        axes(handles.axes4)
        cla reset
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
        if handles.togglebutton7.Value == 0
            plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
            legend('Magn')
        else
            plot(ppms, handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
            legend('Phs. Corr.')
        end
        ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
        end

    end
end



% --- Executes during object creation, after setting all properties.
function CSC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LineBroadFact_Callback(hObject, eventdata, handles)
% hObject    handle to LineBroadFact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LineBroadFact as text
%        str2double(get(hObject,'String')) returns contents of LineBroadFact as a double

[ddim, lbf, tmpdat, tmpdatP, ephci, mepi, FIDdat, imageObj, imageObj2, ...
        phma, dsz, ph1, pivotppm, pivot, imageObj3, siss, pap, FDat3] = datProcLineBroad(handles);

handles = guidata(hObject);
FDat = imageObj.data;
FDat2 = imageObj2.data;
handles.GUIDataAll.FDat = FDat;
handles.GUIDataAll.FDat2 = FDat2;
handles.GUIDataAll.FDat3 = FDat3;
handles.GUIDataAll.lbf = lbf;
guidata(hObject, handles);

if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData_LB',handles.LineBroadFact.String,'.mat'], ''))
    save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData_LB',handles.LineBroadFact.String,'.mat'], ''),'FIDdat', 'FDat2', 'FDat', 'FDat3')
end

% Find SNR value
[SNR, mva] = compSNR(handles, FDat, FDat3);

handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
guidata(hObject, handles);

if isfield(handles.GUIDataAll, 'x2')
    dimsdat = size(handles.GUIDataAll.FDat);
    if handles.togglebutton7.Value == 0 
        fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    else
        fd = handles.GUIDataAll.FDat3(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    end
    x2 = handles.GUIDataAll.x2;
    y2 = handles.GUIDataAll.y2;
    ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);

    if ~isempty(handles.axes4.Children)
        if handles.togglebutton7.Value == 0 
                set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, "String", 'Magn')
            else
                set(handles.axes4.Children, 'YData', real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String))
                set(handles.axes4.Children, 'XData', ppms)
                set(handles.axes4.Legend, 'String', 'Phs. Corr.')
            end

            set(handles.axes4, 'XLim', [min(ppms) max(ppms)])
            set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
    else
        axes(handles.axes4)
        cla reset
        if handles.togglebutton7.Value == 0 
            plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
            legend('Magn')
        else
            plot(ppms, real(handles.GUIDataAll.FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String), 'g')
            legend('Phs. Corr.')
        end
            
        ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
    end

end






% --- Executes during object creation, after setting all properties.
function LineBroadFact_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LineBroadFact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


LineBroadFact_Callback(hObject, [], handles)


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3

LineBroadFact_Callback(hObject, [], handles)


% --- Executes on button press in PlotFID.
function PlotFID_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



try

xx = handles.GUIDataAll.x2;
yy = handles.GUIDataAll.y2;
FIDdat = handles.GUIDataAll.FIDDat.data;


if isempty(strfind(handles.GUIDataAll.FIDDat.Acqp.ACQ_method,'EPSI'))
    timfid = handles.GUIDataAll.FIDDat.Method.PVM_SpecAcquisitionTime;
else
    str = readlines(join([handles.GUIDataAll.CSIgenpath, '\visu_pars'],''));
    mm = [];
    for i = 1:length(str)
        if contains(str(i), 'VisuAcqScanTime')
            mm = char(str(i));
            break
        end
    end
    tt = mm(strfind(mm,'=')+1:end);

    timfid = str2double(tt)/10000;



end
tims = linspace(0, timfid, length(real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))));


if ~isempty(handles.axes2.Children)

    set(handles.axes2.Children, 'YData', real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))
    set(handles.axes2.Children, 'XData', tims)
    set(handles.axes2, 'XLim', [0 timfid])
    set(handles.axes2, 'YLim', [min(real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))) max(real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))])
    set(handles.axes2, 'XDir', 'normal')
    % set(handles.axes2.XLabel, 'String', ' ')
    
    set(handles.axes3.Children, 'YData', imag(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))
    set(handles.axes3.Children, 'XData', tims)
    set(handles.axes3, 'XLim', [0 timfid])
    set(handles.axes3, 'YLim', [min(imag(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))) max(imag(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))))])
    set(handles.axes3, 'XDir', 'normal')
    set(handles.axes3.XLabel, 'String', 'time (ms)')
else

    axes(handles.axes2)  
    cla reset
    plot(tims, real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
    ax = gca;
    legend('Real')
    xlabel('time (ms)')
    
    axes(handles.axes3)
    cla reset
    plot(tims, imag(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
    ax = gca;
    legend('Imag')
    xlabel('time (ms)')
end

catch
end


% --- Executes on button press in togglebutton5.
function togglebutton5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton5

if handles.togglebutton5.Value == 1
    handles = guidata(hObject);
    handles.togglebutton5.BackgroundColor = [0    1    0];
    guidata(hObject, handles); 


    foldpath = handles.GUIDataAll.CSIgenpath;
    if ~isempty(dir(join([foldpath, '/tmp_img/SpectraTimePoint_ColorGrid_',num2str(handles.TimePoints.Value),'_Sat',"*",'.png'],'')))
        
        if handles.radiobutton4.Value == 1
            cg = imread(handles.GUIDataAll.namfil);
        else
            cg = imread(join([foldpath, '/tmp_img/SpectraTimePoint_ColorGrid_',num2str(handles.TimePoints.Value),'_Sat',string(handles.edit8.String),'.png'],''));
        end

        mim = handles.GUIDataAll.CSIVox;

        axes(handles.MainPlot)
        h = get(handles.MainPlot,'Children');
        delete(findobj(h, 'type', 'Patch'));
        delete(findobj(h, 'type', 'Image'));

        handles = guidata(hObject);
        fd = imshow(cg);
        fd.XData = [0 mim(2)];
        fd.YData = [0 mim(3)];
        fd.Parent.XLim = [0,mim(2)];
        fd.Parent.YLim = [0,mim(3)];
        h = get(handles.MainPlot,'Children');
        set(handles.MainPlot,'Children',[h(2:end); h(1)])
        guidata(hObject, handles);

    end


else
    handles = guidata(hObject);
    handles.togglebutton5.BackgroundColor = [0.9400    0.9400    0.9400];
    guidata(hObject, handles);

    handles = guidata(hObject);
    handles.togglebutton1.Value = 0;
    togglebutton1_Callback(hObject, [], handles)

    guidata(hObject, handles);

end


% --- Executes on button press in togglebutton6.
function togglebutton6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton6

if handles.togglebutton4.Value == 1
    handles = guidata(hObject);
    handles.togglebutton4.Value = 0;
    handles.togglebutton4.BackgroundColor = [0.9400    0.9400    0.9400];
%     togglebutton4_Callback(hObject, [], handles);
    guidata(hObject, handles);
end


if handles.togglebutton6.Value == 1
    handles = guidata(hObject);
    handles.togglebutton6.BackgroundColor = [0    1    0];
    guidata(hObject, handles); 

    if isfield(handles, "GUIDataAll") && isfield(handles.GUIDataAll, "FDat")
%         genSpectraAllColorGrid(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String), handles.GUIDataAll.CSIgenpath);
        axes(handles.SpectraPan)
        namfil = genSpectraAllColorGrid_Colormap_MultPeak(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String(handles.TimePoints.Value,:)), handles.GUIDataAll.CSIgenpath, handles.GUIDataAll.ppms, str2num(handles.edit8.String));

%         % THIS IS ONLY FOR PHASE CORRECTION WITH FUMARATE/MALATE
%         namfil = genSpectraAllColorGrid_Colormap_MultPeak_FumarateMaleate(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String(handles.TimePoints.Value,:)), handles.GUIDataAll.CSIgenpath, handles.GUIDataAll.ppms);


        handles.GUIDataAll.namfil = namfil;
        guidata(hObject, handles); 


        im = imread(namfil);
        axes(handles.SpectraPan)

        rati = handles.GUIDataAll.CSISize(2:3)/max(handles.GUIDataAll.CSISize(2:3));

        imm=imshow(im);
        imm.Parent.XLim = [1-0.5   (800+0.5)*rati(1)];
        imm.Parent.YLim = [1-0.5   (800+0.5)*rati(2)];
        imm.XData = [1   800*rati(1)];
        imm.YData = [1   800*rati(2)];


    else
        f = warndlg('You need to load a CSI experiment first if you wanna generate the colored grid! :S','Warning');
        handles = guidata(hObject);
        handles.togglebutton6.Value = 0;
        handles.togglebutton6.BackgroundColor = [0.9400    0.9400    0.9400];
        guidata(hObject, handles);
    end

else
    handles = guidata(hObject);
    handles.togglebutton6.BackgroundColor = [0.9400    0.9400    0.9400];
    guidata(hObject, handles);

    if handles.togglebutton7.Value == 0
        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
    else
        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'_Phased.png'], ''));
    end
    axes(handles.SpectraPan)
    imm=imshow(im);
    imm.Parent.XLim = [1-0.5   800+0.5];
    imm.Parent.YLim = [1-0.5   800+0.5];
    imm.XData = [1   800];
    imm.YData = [1   800];

end


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


togglebutton5_Callback(hObject, [], handles)


% --- Executes on slider movement.
function Phase0Slide_Callback(hObject, eventdata, handles)
% hObject    handle to Phase0Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% disp(handles.Phase0Slide.Value)

if strcmp(handles.togglebutton7.String,'Phase Corr.')


[x2,y2,z,ephci,phma,lbf,dsz,ph1,pivotppm,pivot,tmpdat,siss,ddim,imageObj3,pap,FDat3] = phase0Proc(handles);

handles = guidata(hObject);
handles.GUIDataAll.FDat3 = FDat3;
handles.GUIDataAll.phma = phma;
guidata(hObject, handles);

save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], ''),'FDat3')


% REDUCE TO ONLY TAKE THE ONE VOXEL WE ARE LOOKING AT
% Find SNR value

[SNR, mva] = compSNR(handles, handles.GUIDataAll.FDat, FDat3);

handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
guidata(hObject, handles);


if isfield(handles.GUIDataAll, 'x2')
    dimsdat = size(handles.GUIDataAll.FDat);
    fd = handles.GUIDataAll.FDat3(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    x2 = handles.GUIDataAll.x2;
    y2 = handles.GUIDataAll.y2;
    ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);
    
    axes(handles.axes4)
%     cla reset
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
    
    if handles.togglebutton7.Value == 0 
%         plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
%         legend('Magn')
        set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
        set(handles.axes4.Legend, "String", 'Magn')
    else
%         plot(ppms, real(imageObj3.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String), 'g')
%         legend('Phs. Corr.')
        set(handles.axes4.Children, 'YData', real(FDat3(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String))
        set(handles.axes4.Legend, "String", 'Phs. Corr.')
    end
    
    set(handles.axes4, 'XLim', [min(ppms) max(ppms)])
    set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])

%     ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
%     xlim([min(ppms) max(ppms)])
    ax = gca;
    ax.XDir = 'reverse';
%     set(handles.axes4, 'xlabel', 'ppm')
%     xlabel('ppm')
    
end

handles = guidata(hObject);
handles.Phase0Slide.Value = 0;
handles.GUIDataAll.FDat3 = FDat3;
guidata(hObject, handles);

else
    handles = guidata(hObject);
    handles.Phase0Slide.Value = 0;
    guidata(hObject, handles);
end







% --- Executes during object creation, after setting all properties.
function Phase0Slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phase0Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in togglebutton7.
function togglebutton7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton7


if handles.togglebutton7.Value == 1
    handles = guidata(hObject);
    handles.togglebutton7.BackgroundColor = [0    1    0];
    handles.togglebutton7.String = 'Phase Corr.';
    guidata(hObject, handles); 

else

    handles = guidata(hObject);
    handles.togglebutton7.BackgroundColor = [.9 .9 .9];
    handles.togglebutton7.String = 'Magnitude';
    guidata(hObject, handles); 

end

[SNR, mva] = compSNR(handles, handles.GUIDataAll.FDat, handles.GUIDataAll.FDat3);
    
handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
guidata(hObject, handles);



if handles.togglebutton7.Value == 0
    im = imread(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
else
    im = imread(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'_Phased.png'], ''));
end
axes(handles.SpectraPan)
imm = imshow(im);
imm.Parent.XLim = [1-0.5   800+0.5];
imm.Parent.YLim = [1-0.5   800+0.5];
imm.XData = [1   800];
imm.YData = [1   800];

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function PivotEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PivotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PivotEdit as text
%        str2double(get(hObject,'String')) returns contents of PivotEdit as a double




% --- Executes during object creation, after setting all properties.
function PivotEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PivotEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Phase1Slide_Callback(hObject, eventdata, handles)
% hObject    handle to Phase1Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% disp(handles.Phase1Slide.Value)

if strcmp(handles.togglebutton7.String,'Phase Corr.')

[x2,y2,z,ephci,phma,lbf,dsz,ph1,pivotppm,pivot,tmpdat,siss,ddim,imageObj3,pap,FDat3] = phase1Proc(handles);

handles = guidata(hObject);
handles.GUIDataAll.ph1 = ph1;
handles.GUIDataAll.pivotppm = pivotppm;
handles.GUIDataAll.pivot = pivot;
handles.GUIDataAll.FDat3 = FDat3;
guidata(hObject, handles);

save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhasedData.mat'], ''),'FDat3')

% REDUCE TO ONLY TAKE THE ONE VOXEL WE ARE LOOKING AT
% Find SNR value

[SNR, mva] = compSNR(handles, handles.GUIDataAll.FDat, FDat3);

handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
guidata(hObject, handles);


if isfield(handles.GUIDataAll, 'x2')
    dimsdat = size(handles.GUIDataAll.FDat);
    fd = handles.GUIDataAll.FDat3(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    x2 = handles.GUIDataAll.x2;
    y2 = handles.GUIDataAll.y2;
    ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);

    axes(handles.axes4)
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
%     hold on
    if handles.togglebutton7.Value == 0 
%         plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
%         legend('Magn')
        set(handles.axes4.Children, 'YData', handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String))
        set(handles.axes4.Legend, "String", 'Magn')
    else
%         plot(ppms, real(imageObj3.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String), 'g')
%         legend('Phs. Corr.')
        set(handles.axes4.Children, 'YData', real(imageObj3.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))))*str2double(handles.SFac.String))
        set(handles.axes4.Legend, "String", 'Phs. Corr.')
    end
%     plot([pivotppm(x2,y2,z), pivotppm(x2,y2,z)],[min(real(fd(:)))-max(real(fd(:)))*0.05, max(real(fd(:)))+max(real(fd(:)))*0.05], 'r')
%     hold off
    set(handles.axes4, 'XLim', [min(ppms) max(ppms)])
    set(handles.axes4, 'YLim', [min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
%     ylim([min(real(fd(:)))-max(real(fd(:)))*0.05 max(real(fd(:)))+max(real(fd(:)))*0.05])
%     xlim([min(ppms) max(ppms)])
    
%     ax = gca;
%     ax.XDir = 'reverse';

%     set(handles.axes4, 'xlabel', 'ppm')
%     xlabel('ppm')

end



handles = guidata(hObject);
handles.Phase1Slide.Value = 0;
handles.GUIDataAll.FDat3 = FDat3;
guidata(hObject, handles);
else

    handles = guidata(hObject);
    handles.Phase1Slide.Value = 0;
    guidata(hObject, handles);
end






% --- Executes during object creation, after setting all properties.
function Phase1Slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phase1Slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double

if str2num(handles.edit8.String) < 0
    handles.edit8.String = '0';
elseif str2num(handles.edit8.String) > 100
    handles.edit8.String = '100';
elseif isempty(str2num(handles.edit8.String))
    handles.edit8.String = '0';
end



% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ABL.
function ABL_Callback(hObject, eventdata, handles)
% hObject    handle to ABL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



Phase1Slide_Callback(hObject, [], handles);
