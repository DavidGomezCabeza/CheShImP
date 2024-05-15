function [impix, FOVSizeCSI, FOVPosCSI, FOVPosIm, foldpath, turbo, posIM, posCS, PosGen, xb, yb] = getFOV(handles)

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

    if isempty(strfind(handles.GUIDataAll.FIDDat.Acqp.ACQ_method,'EPSI'))
        xb = ([0 FOVSizeCSI(2) FOVSizeCSI(2) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
        yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
    else % to process EPSI
        xb = ([0 FOVSizeCSI(1) FOVSizeCSI(1) 0]+PosGen(1))*impix(1)/handles.GUIDataAll.MRIImageSize(1);
        yb = ([0 0 FOVSizeCSI(3) FOVSizeCSI(3)]+PosGen(2))*impix(2)/handles.GUIDataAll.MRIImageSize(2);
    end
        


end



