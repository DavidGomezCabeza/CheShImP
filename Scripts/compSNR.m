%% Function to Compute the SNR of each voxel and the maximum

function [SNR, mva] = compSNR(handles, FDat, FDat3)

    perN = 0.1;


    try
        if handles.togglebutton7.Value == 0
            sis = size(FDat);
            FDatFrame = FDat(:,:,:,:,:,:,1);
        else
            sis = size(FDat3);
            FDatFrame = real(FDat3(:,:,:,:,:,:,1));
        end
        mva = max(FDatFrame(:));
        SNR = zeros(1,sis(2)*sis(3));
        SNRAll = zeros(sis(2),sis(3));
        cnt = 1;
        for i = 1:sis(2)
            for j = 1:sis(3)
                SNR(cnt) = max(FDatFrame(:,i,j))/std([FDatFrame(1:round(sis(1)*perN),i,j); FDatFrame(end-round(sis(1)*perN)+1:end,i,j)]);
                cnt = cnt +1;

                SNRAll(i,j) = max(FDatFrame(:,i,j))/std([FDatFrame(1:round(sis(1)*perN),i,j); FDatFrame(end-round(sis(1)*perN)+1:end,i,j)]);
%                 if ismember(mva, FDat(:,i,j))
%                     mvi = [i,j];
%                 end
            end
        end
        save(join([handles.GUIDataAll.CSIgenpath, 'tmp_img\SNR.mat'], ''),'SNR', 'SNRAll')
        
    catch
        SNR = [];
    end
end









