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

% Last Modified by GUIDE v2.5 21-Jul-2023 12:41:20

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




    %% Process data to get right structure
    % Create a RawDataObject, importing the test data

    rawObj = RawDataObject(foldpath, ['fid_proc']);

    
    numSlices = rawObj.Acqp.NI;
    numReps = rawObj.Acqp.NR;
    try
        acqSizes=bruker_getAcqSizes(rawObj.Acqp);
    catch
        acqSizes=bruker_getAcqSizes(rawObj.Acqp, rawObj.Method);
    end



%     fileID = fopen(join([foldpath,'\rawdata.job0'],''),'r');
%     if strcmp(rawObj.Acqp.ACQ_word_size, '_32_BIT')
%         fidFile=fread(fileID, 'int32');
%     elseif strcmp(rawObj.Acqp.ACQ_word_size, '_64_BIT')
%         fidFile=fread(fileID, 'int64');
%     elseif strcmp(rawObj.Acqp.ACQ_word_size, '_16_BIT')
%         fidFile=fread(fileID, 'int16');
%     elseif strcmp(rawObj.Acqp.ACQ_word_size, '_8_BIT')
%         fidFile=fread(fileID, 'int8');
%     end
%     fclose(fileID);
    
try
    fileID = fopen(join([foldpath,'\fid_proc.64'],''),'r');
    fidFile=fread(fileID, 'float64');
    fclose(fileID);

        
    fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));


    % Under TEST -- DAVID
    % Coronal when 3 are taken
    if length(rawObj.Method.PVM_SPackArrSliceOrient)==3
        fidFile2 = fidFile2((acqSizes(1)/2)*acqSizes(2)*acqSizes(3)*numReps*2+1:end);
    end


%     fidFile3 = reshape(fidFile2(1:acqSizes(1)*acqSizes(2)*acqSizes(3)/2), 1, acqSizes(1)/2, acqSizes(2)*acqSizes(3));
    fidFile3 = reshape(fidFile2, 1, acqSizes(1)/2, acqSizes(2)*acqSizes(3)*numReps);

    if length(rawObj.Method.PVM_SPackArrSliceOrient) >= 3 
        rawObj.data{1} = fidFile3;
    end

catch
end


    
%     if acqSizes(2) ~= acqSizes(3)
%         f = warndlg('Carefull! You have a NON squared grid and this is not yet supported! This will generate an error down the line','Warning');
%     end
    
    numPhases=acqSizes(2);
    % Create a FrameDataObject from the imported RawDataObject
    frameObj = FrameDataObject(rawObj);
    
    
    % Under TEST -- DAVID
    % Access coronal for when  3D are acquired
    if length(rawObj.Method.PVM_SPackArrSliceOrient)==3
        frameObj.data(:,:,1,1,:) = fidFile3(1,:,:);
    end
    
    % Create a Cartesian k-space data object from the imported FrameDataObject
    % (Currently this works only with FLASH, MSME, RARE and FISP)
    kdataObj = CKDataObject(frameObj);
%     FDat = kdataObj.data;
    
  
%     sp0 = frameObj.Acqp.ACQ_spatial_phase_0;
%     sp1 = frameObj.Acqp.ACQ_spatial_phase_1;
%     
%     
%     sp0n = 1+((acqSizes(2)-1)/(max(sp0)-min(sp0)))*(sp0 - min(sp0));
%     sp1n = 1+((acqSizes(3)-1)/(max(sp1)-min(sp1)))*(sp1 - min(sp1));
%     
%     
%     siob = size(kdataObj.data);
% 
%     datTmp = zeros(size(kdataObj.data));
%     
%     cnt = 1;
%     for i = 1:acqSizes(2)*acqSizes(3)
%             disp([sp0n(i), sp1n(i)])
%             datTmp(:, sp0n(i), sp1n(i), 1:end, 1:end, 1:end, 1:end, 1:end) = frameObj.data(:, cnt , 1:end, 1:end, 1:end, 1:end, 1:end);
%             cnt = cnt+1;
%             
%     end
% % 
% % frameObj.data = frameObj.data(:, imageObj.Reco.RecoSortMaps(1:end-1)+1, 1:end, 1:end, 1:end, 1:end, 1:end);
% % 
% % 
%     kdataObj.data = datTmp;
% 
% 
    dat = kdataObj.data;
    datsiz = size(dat);
    
    % Here for the data fourier transformed
%     kdataObj2 = CKDataObject(frameObj);
    
%     kdataObj3 = CKDataObject(frameObj);
%     kdataObj3 = kdataObj3.readReco;
    % CHECK THE START OF THE K SPACE, -1 or 0!!!
%     frameObj.data = frameObj.data(:,kdataObj3.Reco.RecoSortMaps(1:end-1)+1);
    
%     mm = zeros(1,16*32);
%     cnt = 1;
%     for i = 1:32
%         mm(1,cnt:cnt+15) = i:32:16*32;
%         cnt = cnt+16;
%     end
% 
%     frameObj.data = frameObj.data(:,mm);


    kdataObj2 = CKDataObject(frameObj);
        
%     kdataObj2.data = datTmp;

    kdataObj2 = kdataObj2.readReco;
    FIDdat = kdataObj2.data;
    % Magnitude
    imageObj=kdataObj2.reco('quadrature');
    imageObj=imageObj.reco('phase_rotate');
    imageObj=imageObj.reco('zero_filling');

    
%     sizdim = size(imageObj.data);
%     imageObj.data = reshape(imageObj.data, sizdim(1), sizdim(2)*sizdim(3));
%     imageObj.data = imageObj.data(:,imageObj.Reco.RecoSortMaps(1:end-1)+1);
%     imageObj.data = imageObj.data(:,255-imageObj.Reco.RecoSortMaps(1:end-1)+1);
%     imageObj.data = reshape(imageObj.data, sizdim(1), sizdim(2),sizdim(3));


%     imageObj=imageObj.reco('FT');
    imageObj.data = fft(imageObj.data);
    imageObj=imageObj.reco('phase_corr_pi');
    imageObj=imageObj.reco('cutoff');
    imageObj=imageObj.reco('scale_phase_channels');
    imageObj=imageObj.reco('sumOfSquares');
    imageObj=imageObj.reco('transposition');
    imageObj.data = flip(imageObj.data,1);

    % Without Magnitude
    imageObj2=kdataObj2.reco('quadrature');
    imageObj2=imageObj2.reco('phase_rotate');
    imageObj2=imageObj2.reco('zero_filling');
%         imageObj2=imageObj2.reco('FT');
    imageObj2.data = fft(imageObj2.data);
    imageObj2=imageObj2.reco('phase_corr_pi');
    imageObj2=imageObj2.reco('cutoff');
    imageObj2=imageObj2.reco('scale_phase_channels');
%     imageObj2=imageObj2.reco('sumOfSquares');
    imageObj2=imageObj2.reco('transposition');
    imageObj2.data = flip(imageObj2.data,1);

    handles = guidata(hObject);
    FDat = imageObj.data;
    FDat2 = imageObj2.data;
    handles.GUIDataAll.FDat = FDat;
    handles.GUIDataAll.FDat2 = FDat2;
    handles.GUIDataAll.Method = rawObj.Method;
    handles.GUIDataAll.numReps = numReps;
    handles.GUIDataAll.FIDDat = kdataObj2;
    guidata(hObject, handles);
 

    if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData.mat'], ''))
        save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData.mat'], ''),'FIDdat', 'FDat2', 'FDat')
    end


    % Find SNR value
    try
        sis = size(FDat);
        FDatFrame = FDat(:,:,:,:,:,:,1);
        mva = max(FDatFrame(:));
        SNR = zeros(1,sis(2)*sis(3));
        cnt = 1;
        for i = 1:sis(2)
            for j = 1:sis(3)
                SNR(cnt) = max(FDatFrame(:,i,j))/std([FDatFrame(1:round(sis(1)*0.1),i,j); FDatFrame(end-round(sis(1)*0.1)+1:end,i,j)]);
                cnt = cnt +1;
%                 if ismember(mva, FDat(:,i,j))
%                     mvi = [i,j];
%                 end
            end
        end

        handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
        guidata(hObject, handles);
        catch
    end


    if ~isfile(join([foldpath, 'tmp_img\CSIProcessedData.mat'], ''))
        Magnitude = FDat;
        Complex = imageObj2.data;
        save(join([foldpath, 'tmp_img\CSIProcessedData.mat'], ''), "Magnitude", "Complex");
    end
    

    if ~isfolder(join([foldpath, 'tmp_img\ProcessedDataCSV'], ''))
        mkdir(join([foldpath, 'tmp_img\ProcessedDataCSV'], ''))
        datsiz2 = size(FDat);
        for i = 1:datsiz2(2)
            for j = 1:datsiz2(3)
                for k = 1:numReps
                    T = array2table([real(imageObj2.data(:,i,j,1,1,1,k)), imag(imageObj2.data(:,i,j,1,1,1,k)), FDat(:,i,j,1,1,1,k)]);
                    T.Properties.VariableNames(1:3) = {'Real','Imaginary','Magnitude'};
                    writetable(T,join([foldpath, 'tmp_img\ProcessedDataCSV\CSIData_X',num2str(i),'_Y',num2str(j),'_T',num2str(k),'.csv'], ''))
                end
            end
        end
    end

    



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


% Time point selection

% if size(datsiz)>3
% %     TPs = datsiz(end);
%     TPs = numReps;
% else
%     TPs = 1;
% end

TPs = numReps;

handles = guidata(hObject);
handles.TimePoints.String = [1:TPs];
guidata(hObject, handles);

    
genSpectraAll(FDat, handles.TimePoints.Value, foldpath)
im = imread(join([foldpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
axes(handles.SpectraPan)
imm=imshow(im);
imm.Parent.XLim = [1-0.5   800+0.5];
imm.Parent.YLim = [1-0.5   800+0.5];
imm.XData = [1   800];
imm.YData = [1   800];

% get the visu parameters for the grid

    try   
        dimIm = [];
        fid = fopen(join([foldpath, 'visu_pars'],''));
        tline = fgetl(fid);
        lineCounter = 1;
        while ischar(tline)
            if contains(tline, 'VisuCoreExtent=', 'IgnoreCase', true)
                tline = fgetl(fid);
                dimIm = tline;
                break;
            end
            % Read next line
            tline = fgetl(fid);
            lineCounter = lineCounter + 1;
        end
        fclose(fid);

        FOVSizeCSI = str2num(dimIm);

        dimIm = [];
        fid = fopen(join([foldpath, 'visu_pars'],''));
        tline = fgetl(fid);
        lineCounter = 1;
        while ischar(tline)
            if contains(tline, '$VisuCorePosition=', 'IgnoreCase', true)
                tline = fgetl(fid);
                dimIm = tline;
                break;
            end
            % Read next line
            tline = fgetl(fid);
            lineCounter = lineCounter + 1;
        end
        fclose(fid);

        FOVPosCSI = str2num(dimIm);

    %     xb = ([0 20 20 0]+19.6797507758513)*256/60;
    %     yb = ([0 0 40 40]+(60-40)-9.93554319775904)*256/60;

        %% THIS IS UNDER TESTING -- DAVID
        

        if handles.FullImage.XLim(2)>5
            
            switch handles.GUIDataAll.Method.PVM_SPackArrSliceOrient
                case 'axial'
                    posIM = abs([handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(2)]);
                    posCS = abs([FOVPosCSI(1) FOVPosCSI(2)]);
                    posCSor = [FOVPosCSI(1) FOVPosCSI(2)];
                    posIMor = [handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(2)];
                case 'sagittal'
                    posIM = abs([handles.GUIDataAll.MRIImagePos(2) handles.GUIDataAll.MRIImagePos(3)]);
                    posCS = abs([FOVPosCSI(2) FOVPosCSI(3)]);
                    posCSor = [FOVPosCSI(2) FOVPosCSI(3)];
                    posIMor = [handles.GUIDataAll.MRIImagePos(2) handles.GUIDataAll.MRIImagePos(3)];
                case 'coronal'
                    posIM = abs([handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(3)]);
                    posCS = abs([FOVPosCSI(1) FOVPosCSI(3)]);
                    posCSor = [FOVPosCSI(1) FOVPosCSI(3)];
                    posIMor = [handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(3)];
            end
            
            PosGen = posIM-posCS;
            impix = size(handles.GUIDataAll.MRIImage);
%             PosGen = handles.GUIDataAll.PosGen;

            if impix(1) == impix(2)
                xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
                yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
            else
                yb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
                xb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
            end
            
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


            

%             PosGen = abs(FOVPosCSI - handles.GUIDataAll.MRIImagePos);
            

            %% David -- THIS IS NOT RIGHT!!!
%             xb = ([0 FOVSizeCSI(3) FOVSizeCSI(3) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);
% 
% %             xb = (([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2));
%             yb = (([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.GUIDataAll.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.GUIDataAll.MRIImageSize(2));

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



        %     figure, 
            mim = size(FDat);

            axes(handles.MainPlot)
%             figure
            fd = imshow(mat2gray(cutMRIIm));
            fd.XData = [0 mim(2)];
            fd.YData = [0 mim(3)];
            fd.Parent.XLim = [0,mim(2)];
            fd.Parent.YLim = [0,mim(3)];
            h = get(handles.MainPlot,'Children');
            set(handles.MainPlot,'Children',[h(2:end); h(1)])
        %     xticks([0:1:mim(2)])
        %     yticks([0:1:mim(3)])
            end
        end
        
    catch

        f = warndlg('No visu_pars file in the same folder than the fid file! Please, fix this if you want to overlay CSI and MRI Image!','Warning');

    end

% Compute the ppm axix for the plots
bw = imageObj2.Method.PVM_SpecSW(1);
bwc = imageObj2.Method.PVM_FrqWorkPpm(1);

   
ppms = flip(linspace(bwc-bw/2, bwc+bw/2, datsiz(1)));
handles.GUIDataAll.ppms = ppms;
guidata(hObject, handles);


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
%     if handles.ClearGUI.Value ==1
%         break;
%     end

    [x,y, but] = ginput(1);
    
    x2=floor(x)+1;
    y2=floor(y)+1;


    if but == 3
        handles.PPM.String = x;
%         if sum(round(ppms) == round(x)) ~= 0
% 
%         else
            [minValue,closestIndex] = min(abs(ppms-x));
            handles.Intens.String = handles.GUIDataAll.currDat(closestIndex);
%         end
%         handles.Intens.String = y;
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
        handles.XYPos.String = join(['X: ?, Y: ?'], '');
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
    fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,handles.TimePoints.Value);

    

    dimsdat = size(imageObj.data);

    try % Check adding the time points, but this might be it!!!
        axes(handles.axes2)  
%         plot(ppms, real(imageObj2.data(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)), 'b')
%         plot(ppms, real(imageObj2.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
        plot(ppms, real(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
        ylim([min(rd(:)) max(rd(:))])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
        legend('Real')
        
        axes(handles.axes3)
%         plot(ppms,imag(imageObj2.data(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)), 'r')
%         plot(ppms,imag(imageObj2.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
        plot(ppms,imag(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
        ylim([min(id(:)) max(id(:))])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
        legend('Imag')
        
        if isnan(str2double(handles.SFac.String))
            axes(handles.axes4)
%             plot(ppms,imageObj.data(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value), 'g')
%             plot(ppms+str2double(handles.CSC.String),imageObj.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))), 'g')
            plot(ppms+str2double(handles.CSC.String),handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:))), 'g')
            ylim([min(fd(:)) max(fd(:))])
            xlim([min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')
            legend('Magn')
            handles = guidata(hObject);
        
%             handles.GUIDataAll.currDat = imageObj.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            handles.GUIDataAll.currDat = handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            guidata(hObject, handles);
        else
            axes(handles.axes4)
%             plot(ppms,imageObj.data(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
%             plot(ppms+str2double(handles.CSC.String),imageObj.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
            plot(ppms+str2double(handles.CSC.String),handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
            ylim([min(fd(:)) max(fd(:))])
            xlim([min(ppms+str2double(handles.CSC.String)) max(ppms+str2double(handles.CSC.String))])
            ax = gca;
            ax.XDir = 'reverse';
            xlabel('ppm')
            legend('Magn')
            handles = guidata(hObject);
       
%             handles.GUIDataAll.currDat = imageObj.data(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
            handles.GUIDataAll.currDat = handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
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
    im = imread(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
    axes(handles.SpectraPan)
    imshow(im)

    % Find SNR value
    try
        FDat = handles.GUIDataAll.FDat;
        sis = size(FDat);
        FDatFrame = FDat(:,:,:,:,:,:,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        mva = max(FDatFrame);
        SNR = zeros(1,sis(2)*sis(3));
        cnt = 1;
        for i = 1:sis(2)
            for j = 1:sis(3)
                SNR(cnt) = max(FDatFrame(:,i,j))/std([FDatFrame(1:round(sis(1)*0.1),i,j); FDatFrame(end-round(sis(1)*0.1)+1:end,i,j)]);
                cnt = cnt +1;
%                 if ismember(mva, FDat(:,i,j))
%                     mvi = [i,j];
%                 end
            end
        end

        handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
        guidata(hObject, handles);
        catch
    end

    handles = guidata(hObject);
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
        f = warndlg('DICOM file NOT selected! Image will be reconstructed from FID file. This is NOT the preferred option! If more than one plane is acquired, only the last one acquired will be displayed!');
        waitfor(f)

        filess = dir(fullfile(join([foldpath, '\pdata\1'],'')));
        for k = 3:size(filess)
            if ~isfile(join([foldpath, filess(k).name],''))
                copyfile(join([foldpath, '\pdata\1\', filess(k).name],''), join([foldpath, '\', filess(k).name],''))
            end
        end

        kdataObj = CKDataObjectIMAGE(foldpath);
        kdataObj = kdataObj.readReco;
        imageObj=kdataObj.reco('all', 'image');

        handles = guidata(hObject);
        handles.edit2.String = join([foldpath, fidfile],'');
        Y = mat2gray(imageObj.data(:,:,1,1,1,end)');
        axes(handles.FullImage)
        imshow(Y,[])
       
        handles.GUIDataAll.MRIImage = Y;
        guidata(hObject, handles);

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
        
        if isfile(join([foldpath, 'visu_pars'],''))
            dimIm = [];
            fid = fopen(join([foldpath, 'visu_pars'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, 'VisuCoreExtent=', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);

            FOVSizeIm = str2num(dimIm);

            dimIm = [];
            fid = fopen(join([foldpath, 'visu_pars'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, '$VisuCorePosition=', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);

            FOVPosIm = str2num(dimIm);

            dimIm = [];
            fid = fopen(join([foldpath, 'acqp'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, 'ACQ_protocol_name', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);

            turbo = contains(lower(dimIm), 'turbo');


        else 
            dimIm = [];
            fid = fopen(join([foldpath, '../visu_pars'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, 'VisuCoreExtent=', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);

            FOVSizeIm = str2num(dimIm);

            dimIm = [];
            fid = fopen(join([foldpath, '../visu_pars'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, '$VisuCorePosition=', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);

            FOVPosIm = str2num(dimIm);

            
            dimIm = [];
            fid = fopen(join([foldpath, '../../../acqp'],''));
            tline = fgetl(fid);
            lineCounter = 1;
            while ischar(tline)
            %     disp(tline)
                if contains(tline, '$ACQ_protocol_name', 'IgnoreCase', true)
                    tline = fgetl(fid);
                    dimIm = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fid);
                lineCounter = lineCounter + 1;
            end
            fclose(fid);
            
            turbo = contains(lower(dimIm), 'turbo');

        end
        
        
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
        impix = size(handles.GUIDataAll.MRIImage);
        FOVSizeCSI = handles.GUIDataAll.CSISize;
        FOVPosCSI = handles.GUIDataAll.CSIPos;
        
        


        switch handles.GUIDataAll.Method.PVM_SPackArrSliceOrient
            case 'axial'
                posIM = abs([handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(2)]);
                posCS = abs([FOVPosCSI(1) FOVPosCSI(2)]);
                posCSor = [FOVPosCSI(1) FOVPosCSI(2)];
                posIMor = [handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(2)];
            case 'sagittal'
                posIM = abs([handles.GUIDataAll.MRIImagePos(2) handles.GUIDataAll.MRIImagePos(3)]);
                posCS = abs([FOVPosCSI(2) FOVPosCSI(3)]);
                posCSor = [FOVPosCSI(2) FOVPosCSI(3)];
                posIMor = [handles.GUIDataAll.MRIImagePos(2) handles.GUIDataAll.MRIImagePos(3)];
            case 'coronal'
                posIM = abs([handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(3)]);
                posCS = abs([FOVPosCSI(1) FOVPosCSI(3)]);
                posCSor = [FOVPosCSI(1) FOVPosCSI(3)];
                posIMor = [handles.GUIDataAll.MRIImagePos(1) handles.GUIDataAll.MRIImagePos(3)];
        end

        if (posCSor(1)<0 && posIMor(1)>0) || (posCSor(1)>0 && posIMor(1)<0) %turbo == 1
%             xb = impix(1) - xb;

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

        xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
        yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);

        % Check if it is using the turbo localiser so we can flip the x
        % axis

        
        





        %% THIS IS UNDER TESTING -- DAVID
%         PosGen = abs(FOVPosCSI - FOVPosIm);

        %% David -- THIS IS NOT RIGHT!!!
%         xb = abs(([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2));
%         xb = ([0 FOVSizeCSI(3) FOVSizeCSI(3) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);

%         yb = abs(([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.GUIDataAll.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.GUIDataAll.MRIImageSize(2));

        axes(handles.FullImage)
        hold on 
        patch(xb,yb,'g', 'EdgeColor', 'none','FaceAlpha',.4);
        
        
        MRIIm = handles.GUIDataAll.MRIImage;
%         cutMRIIm = MRIIm(min(yb):max(yb), min(xb):max(xb));
        cutMRIIm = MRIIm(max([min(yb) 1]):max(yb), max([min(xb) 1]):max(xb));
        imwrite(mat2gray(cutMRIIm), join([handles.GUIDataAll.CSIgenpath, '\tmp_img\MRICut.png'], '')); 
        handles = guidata(hObject);
        handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
        handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
        guidata(hObject, handles);
    %     figure, 
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
    
%     handles.togglebutton1.Value = 1;
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
%         
% 
%         impix = size(handles.MRIImage);
%         FOVSizeCSI = handles.CSISize;
%         FOVPosCSI = handles.CSIPos;
%         FOVPosIm = handles.MRIImagePos;
        foldpath = handles.GUIDataAll.CSIgenpath;
% 
        imOver = imread(join([foldpath, '/tmp_img/SpectraTimePoint_',num2str(handles.TimePoints.Value),'.png'], ''));
%         
%         %% THIS IS UNDER TESTING -- DAVID
%         PosGen = abs(FOVPosCSI - FOVPosIm);
% 
%         xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.MRIImageSize(2);
%         yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.MRIImageSize(2);
% 
        mim = handles.GUIDataAll.CSIVox;
% 
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
        foldpath = handles.GUIDataAll.CSIgenpath;

        if ~isfile(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'],''))            
            system(['python "',foldpath,'/tmp_img/',num2str(handles.TimePoints.Value),'-MRICSIOverlayGen.py"']);
        end
        
        impix = size(handles.GUIDataAll.MRIImage);
        FOVSizeCSI = handles.GUIDataAll.CSISize;
        FOVPosCSI = handles.GUIDataAll.CSIPos;
        FOVPosIm = handles.GUIDataAll.MRIImagePos;
        foldpath = handles.GUIDataAll.CSIgenpath;

        imOver = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'], ''));
        

        if ~isfile(join([foldpath, '\tmp_img\MRICSIOverlay_FullImage.png'], ''))
        % Generate overlay with full image
            imOver2 = imread(join([foldpath, '/tmp_img/MRICSIOverlay_',num2str(handles.TimePoints.Value),'.png'], ''));
            imwrite(mat2gray(handles.GUIDataAll.MRIImage), join([foldpath, '\tmp_img\MRIFullImage.png'], '')); 
    
            % To make sure you get the grid edges in yellow instead of just
            % black
            imOver2(1,:,1) = 237;
            imOver2(1,:,2) = 177;
            imOver2(1,:,3) = 32;
    
            imOver2(:,1,1) = 237;
            imOver2(:,1,2) = 177;
            imOver2(:,1,3) = 32;
    
            imOver2(:,end,1) = 237;
            imOver2(:,end,2) = 177;
            imOver2(:,end,3) = 32;
    
            imcut = handles.GUIDataAll.MRIImage(round(handles.GUIDataAll.xIndCut), round(handles.GUIDataAll.yIndCut));
            imorg = imresize(cat(3, im2uint16(mat2gray(handles.GUIDataAll.MRIImage)), im2uint16(mat2gray(handles.GUIDataAll.MRIImage)), ...
                    im2uint16(mat2gray(handles.GUIDataAll.MRIImage))), size(handles.GUIDataAll.MRIImage)*8);
    
            imorg(round(handles.GUIDataAll.xIndCut(1))*8-7:round(handles.GUIDataAll.xIndCut(end))*8, ...
                    round(handles.GUIDataAll.yIndCut(1))*8-7:round(handles.GUIDataAll.yIndCut(end))*8,:) = imresize(im2uint16(imOver2), size(imcut)*8);
            
            imwrite(imorg, join([foldpath, '\tmp_img\MRICSIOverlay_FullImage.png'], ''));
        end

        %% THIS IS UNDER TESTING -- DAVID
        PosGen = abs(FOVPosCSI - FOVPosIm);

        xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);
        yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.GUIDataAll.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
% 
%         axes(handles.FullImage)
%         hold on 
%         patch(xb,yb,'g', 'EdgeColor', 'none','FaceAlpha',.4);
        
        
%         cutMRIIm = MRIIm(min(yb):max(yb), min(xb):max(xb));
    %     figure, 
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
        f = warndlg('To overlay the CSI image to the MRI DICOM image you need to load both files! :S','Warning');
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
        
        impix = size(handles.GUIDataAll.MRIImage);
        FOVSizeCSI = handles.GUIDataAll.CSISize;
        FOVPosCSI = handles.GUIDataAll.CSIPos;
        FOVPosIm = handles.GUIDataAll.MRIImagePos;
        foldpath = handles.GUIDataAll.CSIgenpath;
        turbo = handles.GUIDataAll.turbo;

        switch handles.GUIDataAll.Method.PVM_SPackArrSliceOrient
            case 'axial'
                posIM = abs([FOVPosIm(1) FOVPosIm(2)]);
                posCS = abs([FOVPosCSI(1) FOVPosCSI(2)]);
            case 'sagittal'
                posIM = abs([FOVPosIm(2) FOVPosIm(3)]);
                posCS = abs([FOVPosCSI(2) FOVPosCSI(3)]);
            case 'coronal'
                posIM = abs([FOVPosIm(1) FOVPosIm(3)]);
                posCS = abs([FOVPosCSI(1) FOVPosCSI(3)]);
        end

        PosGen = posIM-posCS;

        xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
        yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);

        % Check if it is using the turbo localiser so we can flip the x
        % axis

%         if turbo == 1
% %             xb = impix(1) - xb;
% 
%             Y = flip(dicomread(join([foldpath, fidfile],'')),2);
%             axes(handles.FullImage)
%             imshow(Y,[])
%             
%             handles.GUIDataAll.MRIImage = Y;
%             guidata(hObject, handles);
% 
%         end



        %% THIS IS UNDER TESTING -- DAVID
%         PosGen = abs(FOVPosCSI - FOVPosIm);
% 
%         %% David -- THIS IS NOT RIGHT!!!
%         xb = ([0 FOVSizeCSI(3) FOVSizeCSI(3) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);
% 
% %         xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+abs(PosGen(3)))*impix(1)/handles.GUIDataAll.MRIImageSize(2);
%         yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+(handles.GUIDataAll.MRIImageSize(2)-FOVSizeCSI(3))-abs(PosGen(1)))*impix(2)/handles.GUIDataAll.MRIImageSize(2);

%         axes(handles.FullImage)
%         hold on 
%         patch(xb,yb,'g', 'EdgeColor', 'none','FaceAlpha',.4);
        
        
        MRIIm = handles.GUIDataAll.MRIImage;
%         cutMRIIm = MRIIm(min(yb):max(yb), min(xb):max(xb));
        cutMRIIm = MRIIm(max([min(yb) 1]):max(yb), max([min(xb) 1]):max(xb));
        imwrite(mat2gray(cutMRIIm), join([foldpath, '\tmp_img\MRICut.png'], '')); 

        handles = guidata(hObject);
        handles.GUIDataAll.xIndCut = max([min(yb) 1]):max(yb);
        handles.GUIDataAll.yIndCut = max([min(xb) 1]):max(xb);
        guidata(hObject, handles);
    %     figure, 
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
        f = warndlg('To overlay the CSI image to the MRI DICOM image you need to load both files! :S','Warning');
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

handles.ClearGUI.Value = 0;
delete(handles.MainPlot.Children);
handles.DataFolderPath.String = '';
handles.edit2.String = '';
handles.togglebutton1.Value = 0;
handles.togglebutton1.BackgroundColor = [0.9400    0.9400    0.9400];
delete(handles.FullImage.Children)
delete(handles.SpectraPan.Children)
delete(handles.axes2.Children)
delete(handles.axes3.Children)
delete(handles.axes4.Children)
guidata(hObject, handles);



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
        fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        x2 = handles.GUIDataAll.x2;
        y2 = handles.GUIDataAll.y2;
        ppms = handles.GUIDataAll.ppms;

        axes(handles.axes4)
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
        plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
        ylim([min(fd(:)) max(fd(:))])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
        legend('Magn')

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
%             handles.ClearGUI.Value = 1;
    
            delete(join([handles.GUIDataAll.CSIgenpath, '\tmp_img\*'],''))
            closereq(); 
%             pause(1)
            
%             handles.ClearGUI.Value = 0;
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
        genSpectraAllColorGrid_Colormap(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String(handles.TimePoints.Value,:)), handles.GUIDataAll.CSIgenpath);

        im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_ColorGrid_',string(handles.TimePoints.Value),'.png'], ''));
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

    im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
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
        fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
        x2 = handles.GUIDataAll.x2;
        y2 = handles.GUIDataAll.y2;
        ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);


        axes(handles.axes4)
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
        plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
        ylim([min(fd(:)) max(fd(:))])
        xlim([min(ppms) max(ppms)])
        ax = gca;
        ax.XDir = 'reverse';
        xlabel('ppm')
        legend('Magn')

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

ddim = size(handles.GUIDataAll.FDat);
% lbf = exp(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1)));
if handles.radiobutton2.Value == 0
%     lbf = exp(linspace(0,str2num(handles.LineBroadFact.String),ddim(1)));
    lbf = exp(flip(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1))));
else
    lbf = exp(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1)));
end

tmpdat = handles.GUIDataAll.FIDDat;

mepi = zeros(1,ddim(2)*ddim(3)*ddim(end));
cnt = 1;
if handles.radiobutton3.Value == 1
    if ddim(2) == ddim(3)
        for i = 1:ddim(2)
            for j = 1:ddim(3)
                if length(ddim)>3
                    for k = 1:ddim(end)
                        mepi(cnt) = find(max(tmpdat.data(:,i,j,:,:,:,k)) == tmpdat.data(:,i,j,:,:,:,k));
                        cnt = cnt+1;
                    end
                else
                    mepi(cnt) = find(max(tmpdat.data(:,i,j)) == tmpdat.data(:,i,j));
                    cnt = cnt+1;
                end
            end
        end
    else
        for i = 1:ddim(2)
            for j = 1:ddim(3)
                if length(ddim)>3
                    for k = 1:ddim(end)
                        mepi(cnt) = find(max(tmpdat.data(:,j,i,:,:,:,k)) == tmpdat.data(:,j,i,:,:,:,k));
                        cnt = cnt+1;
                    end
                else
                    mepi(cnt) = find(max(tmpdat.data(:,j,i)) == tmpdat.data(:,j,i));
                    cnt = cnt+1;
                end
            end
        end
    end

    tmp1 = exp(flip(-linspace(0,str2num(handles.LineBroadFact.String),round(mean(mepi)))));
    tmp2 = exp(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1)-round(mean(mepi))+1));
    lbf = [tmp1, tmp2(2:end)];

end


if ddim(2) == ddim(3)
    for i = 1:ddim(2)
        for j = 1:ddim(3)
            if length(ddim)>3
                for k = 1:ddim(end)
                    tmpdat.data(:,i,j,:,:,:,k) = tmpdat.data(:,i,j,:,:,:,k).*lbf';
                end
            else
                tmpdat.data(:,i,j) = tmpdat.data(:,i,j).*lbf';
            end
        end
    end
else
    for i = 1:ddim(2)
        for j = 1:ddim(3)
            if length(ddim)>3
                for k = 1:ddim(end)
                    tmpdat.data(:,j,i,:,:,:,k) = tmpdat.data(:,j,i,:,:,:,k).*lbf';
                end
            else
                tmpdat.data(:,j,i) = tmpdat.data(:,j,i).*lbf';
            end
        end
    end
end

FIDdat = tmpdat.data;

imageObj=tmpdat.reco('quadrature');
imageObj=imageObj.reco('phase_rotate');
imageObj=imageObj.reco('zero_filling');
imageObj.data = fft(imageObj.data);
imageObj=imageObj.reco('phase_corr_pi');
imageObj=imageObj.reco('cutoff');
imageObj=imageObj.reco('scale_phase_channels');
imageObj=imageObj.reco('sumOfSquares');
imageObj=imageObj.reco('transposition');
imageObj.data = flip(imageObj.data,1);

imageObj2=tmpdat.reco('quadrature');
imageObj2=imageObj2.reco('phase_rotate');
imageObj2=imageObj2.reco('zero_filling');
imageObj2.data = fft(imageObj2.data);
imageObj2=imageObj2.reco('phase_corr_pi');
imageObj2=imageObj2.reco('cutoff');
imageObj2=imageObj2.reco('scale_phase_channels');
imageObj2=imageObj2.reco('transposition');
imageObj2.data = flip(imageObj2.data,1);


handles = guidata(hObject);
FDat = imageObj.data;
FDat2 = imageObj2.data;
handles.GUIDataAll.FDat = FDat;
handles.GUIDataAll.FDat2 = FDat2;
guidata(hObject, handles);

if ~isfile(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData_LB',handles.LineBroadFact.String,'.mat'], ''))
    save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\AllData_LB',handles.LineBroadFact.String,'.mat'], ''),'FIDdat', 'FDat2', 'FDat')
end

% Find SNR value
try
    sis = size(FDat);
    FDatFrame = FDat(:,:,:,:,:,:,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    mva = max(FDatFrame(:));
    SNR = zeros(1,sis(2)*sis(3));
    cnt = 1;
    for i = 1:sis(2)
        for j = 1:sis(3)
            SNR(cnt) = max(FDatFrame(:,i,j))/std([FDatFrame(1:round(sis(1)*0.1),i,j); FDatFrame(end-round(sis(1)*0.1)+1:end,i,j)]);
            cnt = cnt +1;
        end
    end

    handles.text11.String = join(['Max SNR: ', num2str(max(SNR))], '');
    guidata(hObject, handles);
catch
end

if isfield(handles.GUIDataAll, 'x2')
    dimsdat = size(handles.GUIDataAll.FDat);
    fd = handles.GUIDataAll.FDat(:,:,:,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)));
    x2 = handles.GUIDataAll.x2;
    y2 = handles.GUIDataAll.y2;
    ppms = handles.GUIDataAll.ppms+str2double(handles.CSC.String);
    
    rd = real(handles.GUIDataAll.FDat2(:,:,:,1,1,1,handles.TimePoints.Value));
    id = imag(handles.GUIDataAll.FDat2(:,:,:,1,1,1,handles.TimePoints.Value));

    axes(handles.axes2)  
    plot(ppms, real(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
    ylim([min(rd(:)) max(rd(:))])
    xlim([min(ppms) max(ppms)])
    ax = gca;
    ax.XDir = 'reverse';
    xlabel('ppm')
    legend('Real')

    
    axes(handles.axes3)
    plot(ppms,imag(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
    ylim([min(id(:)) max(id(:))])
    xlim([min(ppms) max(ppms)])
    ax = gca;
    ax.XDir = 'reverse';
    xlabel('ppm')
    legend('Imag')


    axes(handles.axes4)
%         plot(ppms, handles.GUIDataAll.FDat(:,dimsdat(2)-x2+1,dimsdat(3)-y2+1,1,1,1,handles.TimePoints.Value)*str2double(handles.SFac.String), 'g')
    plot(ppms, handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
    ylim([min(fd(:)) max(fd(:))])
    xlim([min(ppms) max(ppms)])
    ax = gca;
    ax.XDir = 'reverse';
    xlabel('ppm')
    legend('Magn')

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

axes(handles.axes2)  
plot(real(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
ax = gca;
legend('Real')


axes(handles.axes3)
plot(imag(FIDdat(:,yy,xx,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
ax = gca;
legend('Imag')

% figure, 
% subplot(2,1,1);
% plot(real(FIDdat(:,xx,yy,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
% subplot(2,1,2); 
% plot(imag(FIDdat(:,xx,yy,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')

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
    if isfile(join([foldpath, '/tmp_img/SpectraTimePoint_ColorGrid_',num2str(handles.TimePoints.Value),'.png'],''))

        if handles.radiobutton4.Value == 1
            cg = imread(handles.GUIDataAll.namfil);
        else
            cg = imread(join([foldpath, '/tmp_img/SpectraTimePoint_ColorGrid_',num2str(handles.TimePoints.Value),'.png'],''));
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
        namfil = genSpectraAllColorGrid_Colormap_MultPeak(handles.GUIDataAll.FDat, str2num(handles.TimePoints.String(handles.TimePoints.Value,:)), handles.GUIDataAll.CSIgenpath, handles.GUIDataAll.ppms);

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

    im = imread(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SpectraTimePoint_',string(handles.TimePoints.Value),'.png'], ''));
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
