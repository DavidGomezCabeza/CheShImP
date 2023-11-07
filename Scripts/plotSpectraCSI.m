function plotSpectraCSI(FDat, tf, foldpath)


%     if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.png'], ''))
    
        fsiz = size(FDat);
        totDat = real(FDat(:,:,:,1,1,1,:));
    
        gensiz = size(FDat);
        subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);
    
%         ff=figure('Visible','off');

        ff=figure();
            mm = zeros(1,gensiz(2)*gensiz(3));
            cnt = 1;
            for i = 1:gensiz(2)
                mm(1,cnt:cnt+gensiz(3)-1) = i:gensiz(2):gensiz(2)*gensiz(3);
                cnt = cnt+gensiz(3);
            end
        
        count = 1;
        for i = 1:fsiz(2)    
            for j = 1:fsiz(3)
%                 Double check if this is fliped or not!!!
%                 plot((real((FDat(:,(fsiz(2)+1)-j,(fsiz(2)+1)-i,1,1,1,tf)))), 'g','LineWidth',1)
                if gensiz(2)==gensiz(3) % squared grid
                    subplot(gensiz(2),gensiz(3),count)
                    count= count+1;
                    plot((real((FDat(:,j,i,1,1,1,tf)))), 'g','LineWidth',1)
                else % Rectangular grid
                    subplot(gensiz(3),gensiz(2),mm(count))
                    count= count+1;
                    
                    plot(permute(real((FDat(:,i,j,1,1,1,tf))), [1,3,2,4,5,6,7,8]), 'g','LineWidth',1)
                end
                ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
                xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])
                set(gcf, 'color', 'none');   
                set(gca, 'color', 'none');
                set(gca,'XTick',[], 'YTick', [])
                set(gca, 'XColor','#EDB120', 'YColor','#EDB120')
    
            end
        end
        set(gca,'XTick',[], 'YTick', [])

        
        
    
%         saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.png'], ''))
%         saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.fig'], ''))
        
%         close(ff)
    
%     end

end