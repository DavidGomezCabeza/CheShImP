

function genSpectraAllColorGrid(FDat, tf, foldpath)


    if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''))
    
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
                else % Rectangular grid
                    subplot(gensiz(3),gensiz(2),mm(count))
                    count= count+1;
                    plot(permute(real((FDat(:,i,j,1,1,1,tf))), [1,3,2,4,5,6,7,8]), 'black','LineWidth',1)
                end
                ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
                xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])
%                 set(gcf, 'color', 'none');   

%                 if max(FDat(:,i,j,1,1,1,tf)) >= max(FDat(:))*0.25
%                     set(gca, 'color', 'green');
%                 elseif max(FDat(:,i,j,1,1,1,tf)) <= max(FDat(:))*0.25 && max(FDat(:,i,j,1,1,1,tf)) >= max(FDat(:))*0.1
%                     set(gca, 'color', 'blue');
%                 else
%                     set(gca, 'color', 'red');
%                 end
                updat = max(FDat(:))-mean(FDat(:));
                if max(FDat(:,i,j,1,1,1,tf)) >= mean(FDat(:))+(updat*0.30)
                    set(gca, 'color', '#b5ff84');
                    set(gcf, 'color', '#b5ff84'); 
                elseif max(FDat(:,i,j,1,1,1,tf)) <= mean(FDat(:))+(updat*0.30) && max(FDat(:,i,j,1,1,1,tf)) >= mean(FDat(:))+(updat*0.10)
                    set(gca, 'color', '#ffe18a');
                    set(gcf, 'color', '#ffe18a'); 
                else
                    set(gca, 'color', '#ff837a');
                    set(gcf, 'color', '#ff837a'); 
                end
                
                
                set(gca,'XTick',[], 'YTick', [])
                set(gca, 'XColor','#e0dede', 'YColor','#e0dede')
            end
        end
        set(gca,'XTick',[], 'YTick', [])
        set(gcf, 'InvertHardcopy', 'off');
       
%         daspect(ff.Children, [gensiz(2) gensiz(3) 1])
%         set(gcf, 'Units', 'Inches', 'PaperUnits', 'Inches', 'PaperSize', [700.25, 9.125])

        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.png'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.svg'], ''))
        
        close(ff)
        
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');

    else
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');
    end

end


