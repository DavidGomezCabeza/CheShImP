%% Function to load the FID data and process it to generate spectra in magnitude mode

function [rawObj, numSlices, numReps, acqSizes, fidFile2, fidFile3, numPhases, frameObj, ...
    kdataObj, dat, datsiz, kdataObj2, FIDdat, imageObj, imageObj2] = processMagDat(foldpath)


    rawObj = RawDataObject(foldpath, ['fid_proc']);

    
    numSlices = rawObj.Acqp.NI;
    numReps = rawObj.Acqp.NR;
    try
        acqSizes=bruker_getAcqSizes(rawObj.Acqp);
    catch
        acqSizes=bruker_getAcqSizes(rawObj.Acqp, rawObj.Method);
    end


    if ~isempty(strfind(rawObj.Acqp.ACQ_method,'EPSI'))
        acqSizes = [acqSizes(1)/rawObj.Method.PVM_EncMatrix(1)*rawObj.Method.SpecSize, rawObj.Method.PVM_EncMatrix(1), rawObj.Method.PVM_EncMatrix(2)];
%         [paramStruct.ACQ_size(1)/paramStruct.PVM_EncMatrix(1)*paramStruct.ACQ_size(2), paramStruct.PVM_EncMatrix(1), paramStruct.PVM_EncMatrix(2)];
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
    
 if isempty(strfind(rawObj.Acqp.ACQ_method,'EPSI'))
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
 else

fidFile2 = [];
fidFile3 = [];
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
    
 
    dat = kdataObj.data;
    datsiz = size(dat);
    

    kdataObj2 = CKDataObject(frameObj);
        
%     kdataObj2.data = datTmp;

    kdataObj2 = kdataObj2.readReco;

    if ~isempty(strfind(rawObj.Acqp.ACQ_method,'EPSI')) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS MIGHT NOT BE POSSIBLE
        kdataObj.data = permute(kdataObj.data,[1,3,2]);
        kdataObj2.data = permute(kdataObj2.data,[1,3,2]);
    end


    FIDdat = kdataObj2.data;
    % Magnitude
    imageObj=kdataObj2.reco('quadrature');
    imageObj=imageObj.reco('phase_rotate');
    imageObj=imageObj.reco('zero_filling');

  
%     imageObj=imageObj.reco('FT');
    imageObj.data = fft(imageObj.data); % NEED TO CHECK IF FOR EPSI WE HAVE SPATIAL DIMENSIONS FOURIER TRANSFORMED OR NOT
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

















end










