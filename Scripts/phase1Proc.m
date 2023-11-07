%% Main script to process phase 1 correction for individual spectra

function [x2,y2,z,ephci,phma,lbf,dsz,ph1,pivotppm,pivot,tmpdat,siss,ddim,imageObj3,pap,FDat3] = phase1Proc(handles)


    x2 = handles.GUIDataAll.x2;
    y2 = handles.GUIDataAll.y2;
    z = handles.TimePoints.Value;
    
    ephci = handles.GUIDataAll.ephci;
    
    phma = handles.GUIDataAll.phma;
    
    
    dsz = handles.GUIDataAll.dsz;
    
    ph1 = handles.GUIDataAll.ph1;
    
    ph1(x2,y2,z) = ph1(x2,y2,z) + handles.Phase1Slide.Value;
    
    
    pivotppm = handles.GUIDataAll.pivotppm;
    pivotppm(x2,y2,z) = str2num(handles.PivotEdit.String);
    
    pivot = handles.GUIDataAll.pivot;
    pivot(x2,y2,z) = length(handles.GUIDataAll.ppms((handles.GUIDataAll.ppms <= pivotppm(x2,y2,z))))+1;
    
    try
        lbf = handles.GUIDataAll.lbf;
    catch
        lbf = 1;
    end
    
    save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\PhaseParams.mat'], ''),'ephci', 'phma', 'dsz', 'ph1', 'pivotppm', 'pivot')
    
    tmpdat = handles.GUIDataAll.FIDDat;
    

    % THIS IS FOR PHASE CORRECTION!!! You have to remove phase_corr_pi!
    siss = size(tmpdat.data);
    for i = 1:siss(2)
        for j = 1:siss(3)
            if length(siss) >3
                for z = 1:siss(7)
                    tmpdat.data(:,i,j,1,1,1,z) = [tmpdat.data(ephci+1:end,i,j,1,1,1,z); tmpdat.data(1:ephci,i,j,1,1,1,z)];
                end
            else
                tmpdat.data(:,i,j) = [tmpdat.data(ephci+1:end,i,j); tmpdat.data(1:ephci,i,j)];
            end
        end
    end
    
    
    ddim = size(tmpdat.data);
    
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
        for i = 1:ddim(3)
            for j = 1:ddim(2)
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
    
    
    
    
    imageObj3=tmpdat.reco('quadrature');
    imageObj3=imageObj3.reco('phase_rotate');
    imageObj3=imageObj3.reco('zero_filling');
    
    
    imageObj3.data = fft(imageObj3.data);
    imageObj3=imageObj3.reco('cutoff');
    imageObj3=imageObj3.reco('scale_phase_channels');
    imageObj3=imageObj3.reco('transposition');
    imageObj3.data = flip(imageObj3.data,1);
    
    
    siss = size(imageObj3.data);
    
    if length(siss) > 3
        pap = siss(7);
    else
        pap = 1;
    end
    
    FDat3 = handles.GUIDataAll.FDat3;
    
    z = handles.TimePoints.Value;
    
    
    % imageObj3.data(:,x2,y2,1,1,1,z) = imageObj3.data(:,x2,y2,1,1,1,z) .* exp(sqrt(-1)  *  phma(x2,y2,z) +ph1(x2,y2,z).*((-pivot(x2,y2,z):-pivot(x2,y2,z)+dsz-1)/dsz)'  );
    
    imageObj3.data(:,x2,y2,1,1,1,z) = imageObj3.data(:,x2,y2,1,1,1,z) .* exp( ...
                sqrt(-1) .* (phma(x2,y2,z) + ph1(x2,y2,z) .* ( ...
                    (-pivot(x2,y2,z):-pivot(x2,y2,z)+dsz-1)/dsz)'));
    
    
    
    % Baseline correction from https://es.mathworks.com/matlabcentral/fileexchange/69649-raman-spectrum-baseline-removal
    if handles.ABL.Value == 1
        [~, imageObj3.data(:,x2,y2,1,1,1,z)]=baseline(imageObj3.data(:,x2,y2,1,1,1,z));
    end
    
    
    FDat3(:,x2,y2,1,1,1,z) = imageObj3.data(:,x2,y2,1,1,1,z);



end








