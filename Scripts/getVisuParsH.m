%% Function to extract visualisation parameters from ParaVision 360 for Proton Images

function [FOVSizeIm, FOVPosIm, turbo] = getVisuParsH(foldpath)

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

end