function [ddim, lbf, tmpdat, tmpdatP, ephci, mepi, FIDdat, imageObj, imageObj2, ...
        phma, dsz, ph1, pivotppm, pivot, imageObj3, siss, pap, FDat3] = datProcLineBroad(handles)

    ddim = size(handles.GUIDataAll.FDat);

    if handles.radiobutton2.Value == 0
        lbf = exp(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1)));
    else
        
        lbf = exp(flip(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1))));
    end
    
    tmpdat = handles.GUIDataAll.FIDDat;
    tmpdatP = handles.GUIDataAll.FIDDat;
    ephci = handles.GUIDataAll.ephci;
    
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
    
                            tmpdatP.data(:,i,j,1,1,1,k) = [tmpdatP.data(ephci+1:end,i,j,1,1,1,k); tmpdatP.data(1:ephci,i,j,1,1,1,k)];
                        end
                    else
                        mepi(cnt) = find(max(tmpdat.data(:,i,j)) == tmpdat.data(:,i,j));
                        cnt = cnt+1;
    
                        tmpdatP.data(:,i,j) = [tmpdatP.data(ephci+1:end,i,j); tmpdatP.data(1:ephci,i,j)];
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
    
                            tmpdatP.data(:,j,i,1,1,1,k) = [tmpdatP.data(ephci+1:end,j,i,1,1,1,k); tmpdatP.data(1:ephci,j,i,1,1,1,k)];
                        end
                    else
                        mepi(cnt) = find(max(tmpdat.data(:,j,i)) == tmpdat.data(:,j,i));
                        cnt = cnt+1;
    
                        tmpdatP.data(:,j,i) = [tmpdatP.data(ephci+1:end,j,i); tmpdatP.data(1:ephci,j,i)];
                    end
                end
            end
        end
        tmp1 = exp(flip(-linspace(0,str2num(handles.LineBroadFact.String),round(mean(mepi)))));
        tmp2 = exp(-linspace(0,str2num(handles.LineBroadFact.String),ddim(1)-round(mean(mepi))+1));
        lbf = [tmp1, tmp2(2:end)];
    else
        if ddim(2) == ddim(3)
            for i = 1:ddim(2)
                for j = 1:ddim(3)
                    if length(ddim)>3
                        for k = 1:ddim(end)
                            tmpdatP.data(:,i,j,1,1,1,k) = [tmpdatP.data(ephci+1:end,i,j,1,1,1,k); tmpdatP.data(1:ephci,i,j,1,1,1,k)];
                        end
                    else
                        tmpdatP.data(:,i,j) = [tmpdatP.data(ephci+1:end,i,j); tmpdatP.data(1:ephci,i,j)];
                    end
                end
            end
        else
            for i = 1:ddim(2)
                for j = 1:ddim(3)
                    if length(ddim)>3
                        for k = 1:ddim(end)
                            tmpdatP.data(:,j,i,1,1,1,k) = [tmpdatP.data(ephci+1:end,j,i,1,1,1,k); tmpdatP.data(1:ephci,j,i,1,1,1,k)];
                        end
                    else
                        tmpdatP.data(:,j,i) = [tmpdatP.data(ephci+1:end,j,i); tmpdatP.data(1:ephci,j,i)];
                    end
                end
            end
        end
    
    end
    
    
    if ddim(2) == ddim(3)
        for i = 1:ddim(2)
            for j = 1:ddim(3)
                if length(ddim)>3
                    for k = 1:ddim(end)
                        tmpdat.data(:,i,j,:,:,:,k) = tmpdat.data(:,i,j,:,:,:,k).*lbf';
                        tmpdatP.data(:,i,j,:,:,:,k) = tmpdatP.data(:,i,j,:,:,:,k).*lbf';
                    end
                else
                    tmpdat.data(:,i,j) = tmpdat.data(:,i,j).*lbf';
                    tmpdatP.data(:,i,j) = tmpdatP.data(:,i,j).*lbf';
                end
            end
        end
    else
        for i = 1:ddim(2)
            for j = 1:ddim(3)
                if length(ddim)>3
                    for k = 1:ddim(end)
                        tmpdat.data(:,j,i,:,:,:,k) = tmpdat.data(:,j,i,:,:,:,k).*lbf';
                        tmpdatP.data(:,j,i,:,:,:,k) = tmpdatP.data(:,j,i,:,:,:,k).*lbf';
                    end
                else
                    tmpdat.data(:,j,i) = tmpdat.data(:,j,i).*lbf';
                    tmpdatP.data(:,j,i) = tmpdatP.data(:,j,i).*lbf';
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
    
    
    ephci = handles.GUIDataAll.ephci;
    phma = handles.GUIDataAll.phma;
    dsz = handles.GUIDataAll.dsz;
    ph1 = handles.GUIDataAll.ph1;
    pivotppm = handles.GUIDataAll.pivotppm;
    pivot = handles.GUIDataAll.pivot;
    
    
    imageObj3=tmpdatP.reco('quadrature');
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
    
    for z = 1:pap
        for i = 1:siss(2)
            for j = 1:siss(3)
                    imageObj3.data(:,i,j,1,1,1,z) = imageObj3.data(:,i,j,1,1,1,z) .* exp( ...
                                        sqrt(-1) .* (phma(i,j,z) + ph1(i,j,z) .* ( ...
                                        (-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)'));
                    if ph1(i,j,z) ~= 0
                        % Baseline correction from https://es.mathworks.com/matlabcentral/fileexchange/69649-raman-spectrum-baseline-removal
                        [~, imageObj3.data(:,i,j,1,1,1,z)]=baseline(imageObj3.data(:,i,j,1,1,1,z));
                    end
            end
        end
    end
    FDat3 = imageObj3.data;


end