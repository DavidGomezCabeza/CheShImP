%% Function to extract FOV coordenates for Proton Images

function [impix, FOVSizeCSI, FOVPosCSI, posIM, posCS, posCSor, posIMor] = getFOVCorsH(handles)

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


end