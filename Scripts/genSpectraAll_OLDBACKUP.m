function genSpectraAll(FDat, tf, foldpath)

% FDat

% tf = 1;

if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.png'], ''))

    fsiz = size(FDat);
    % totDat = FDat(:,:,:,1,1,1,tf);
    % totDat = real(FDat(:,:,:,1,1,1,:));
    % totDat = real(fft(FDat(:,:,:,1,1,1,:)));
    
    
    totDat = FDat(:,:,:,1,1,1,:);
    
    
    
%     totDat = real(FDat(:,:,:,1,1,1,:));


    gensiz = size(FDat);
    subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);

%     ff=figure('Visible','off');
    ff=figure();

%     [ha, pos] = tight_subplot(fsiz(2),fsiz(3),0,0);
    % for k = 1:length(ha)
    %     ha(k).Visible = 'off';
    % end
    
    
    count = 1;
    for i = 1:fsiz(2)    
        for j = 1:fsiz(3)
            subplot(gensiz(2),gensiz(3),count)
            count= count+1;
            plot(real((FDat(:,i,j,1,1,1,tf))), 'g','LineWidth',1)

%             axes(ha(count)) 
%             count= count+1;
    %         plot(FDat(:,j,(fsiz(2)+1)-i,1,1,1,tf),'LineWidth',3)
    %         plot(FDat(:,j,i,1,1,1,tf),'LineWidth',3)
    %         plot(real((FDat(:,j,i,1,1,1,tf))),'LineWidth',1)
    %         plot(real((FDat(:,(fsiz(2)+1)-i,(fsiz(2)+1)-j,1,1,1,tf))),'LineWidth',1)
%             plot(real((FDat(:,i,j,1,1,1,tf))), 'g','LineWidth',1)
            
        
            ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
            xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])
            set(gcf, 'color', 'none');   
            set(gca, 'color', 'none');
            set(gca,'XTick',[], 'YTick', [])
            set(gca, 'XColor','#EDB120', 'YColor','#EDB120')

        end
    end
%     set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
    set(gca,'XTick',[], 'YTick', [])
    
    
%     F= getframe(ff);
%     img = F.cdata;
%     imwrite(img, join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.png'], ''));

    
    saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.png'], ''))
    saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_',string(tf),'.fig'], ''))
    
    close(ff)

end

end