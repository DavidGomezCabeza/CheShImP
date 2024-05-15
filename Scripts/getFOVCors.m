%% Function to extract FOV coordenates for CSI and Proton Images

function [posIM, posCS, posCSor, posIMor, xb, yb] = getFOVCors(handles, FOVPosCSI, FOVSizeCSI)


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

    if impix(1) == impix(2)
        if isempty(strfind(handles.GUIDataAll.FIDDat.Acqp.ACQ_method,'EPSI'))
            xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
            yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
        else % This is for EPSI images
            xb = ([0 FOVSizeCSI(1) FOVSizeCSI(1) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
            yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
        end
    else
        if isempty(strfind(handles.GUIDataAll.FIDDat.Acqp.ACQ_method,'EPSI'))
            yb = ([0 FOVSizeCSI(1) FOVSizeCSI(1) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
            xb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
        else % This is for EPSI images

        end
    end

end




