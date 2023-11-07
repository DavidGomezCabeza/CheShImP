%% Function to extract visualisation parameters from ParaVision 360 for CSI

function [FOVSizeCSI, FOVPosCSI] = getVisuPars(handles, foldpath)

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

end










