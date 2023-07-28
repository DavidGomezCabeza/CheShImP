

function genSpectraAllColorGrid_Colormap(FDat, tf, foldpath)


    if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''))
        
        % This assigns a different colour at each pixel almost for sure
%         FDatMaxs = max(FDat);
% %         colmp = colormap(autumn(length(FDatMaxs(:))));
%         colmp = colormap(parula(length(FDatMaxs(:))));
% %         colmp = colormap(turbo(length(FDatMaxs(:))));
% %         colmp = colormap(hot(length(FDatMaxs(:))));
% %         colmp = colormap(cool(length(FDatMaxs(:))));
% %         colmp = colormap(summer(length(FDatMaxs(:))));
% %         colmp = colormap(winter(length(FDatMaxs(:))));
% %         colmp = colormap(spring(length(FDatMaxs(:))));
% %         colmp = colormap(copper(length(FDatMaxs(:))));
%         sortDat = sort(FDatMaxs(:));

        % Groups voxels with comparable maximum peak value
        FDatMaxs = max(FDat);
        FDatMaxsCmp = round(FDatMaxs(:)/min(FDatMaxs(:)));
        colmp = colormap(parula(length(unique(FDatMaxsCmp(:)))));
        colormap(gray)
%         colmp = colormap(autumn(length(unique(FDatMaxsCmp(:)))));
        sortDat = sort(unique(FDatMaxsCmp(:)));




        fsiz = size(FDat);
        totDat = real(FDat(:,:,:,1,1,1,:));
    
        gensiz = size(FDat);
        subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);
    
        ff=figure('Visible','off');
    %     ff=figure();
        
            mm = zeros(1,gensiz(2)*gensiz(3));
            cnt = 1;
            for i = 1:gensiz(2)
                mm(1,cnt:cnt+gensiz(3)-1) = i:gensiz(2):gensiz(2)*gensiz(3);
                cnt = cnt+gensiz(3);
            end
        
        count = 1;
        for i = 1:fsiz(2)    
            for j = 1:fsiz(3)
                if gensiz(2)==gensiz(3) % squared grid
                    subplot(gensiz(2),gensiz(3),count)
                    count= count+1;
                    plot((real((FDat(:,j,i,1,1,1,tf)))), 'black','LineWidth',1)
%                     mxd = max(FDat(:,j,i,1,1,1,tf)); % First case where you have a different colour for each voxel 
                    mxd = round(max(FDat(:,j,i,1,1,1,tf))/min(FDatMaxs(:))); % Second case where you group voxels of comparable intensity
                else % Rectangular grid
                    subplot(gensiz(3),gensiz(2),mm(count))
                    count= count+1;
                    plot(permute(real((FDat(:,i,j,1,1,1,tf))), [1,3,2,4,5,6,7,8]), 'black','LineWidth',1)
%                     mxd = max(FDat(:,i,j,1,1,1,tf)); % First case where you have a different colour for each voxel 
                    mxd = round(max(FDat(:,i,j,1,1,1,tf))/min(FDatMaxs(:))); % Second case where you group voxels of comparable intensity
                end
                ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
                xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])


                set(gca, 'color', colmp(find(mxd == sortDat), :));
                set(gcf, 'color', colmp(find(mxd == sortDat), :)); 

                
                set(gca,'XTick',[], 'YTick', [])
                set(gca, 'XColor','#e0dede', 'YColor','#e0dede')
            end
        end
        set(gca,'XTick',[], 'YTick', [])
        set(gcf, 'InvertHardcopy', 'off');
 

        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.png'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''))
        
        close(ff)
        
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');

    else
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');
    end

end


