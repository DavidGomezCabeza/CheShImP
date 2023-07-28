

FDat2 = reshape(FDat, 128,64,15);

totDat = real(FDat(:,:,:,1,1,1,:));
ff=figure('Visible','off');
[ha, pos] = tight_subplot(fsiz(2),fsiz(3),0,0);



count = 1;
% cr = [1:2:64, 2:2:64];
% cr = [1:3:64, 2:3:64, 3:3:64];
% cr = [1:4:64, 2:4:64, 3:4:64, 4:4:64];
% cr = [1:5:64, 2:5:64, 3:5:64, 4:5:64, 5:5:64];
% cr = [1:6:64, 2:6:64, 3:6:64, 4:6:64, 5:6:64, 6:6:64];
% cr = [1:7:64, 2:7:64, 3:7:64, 4:7:64, 5:7:64, 6:7:64, 7:7:64];
cr = [1:8:64, 2:8:64, 3:8:64, 4:8:64, 5:8:64, 6:8:64, 7:8:64, 8:8:64];


for i = 1:64    
%     for j = 1:fsiz(3)
        axes(ha(count)) 
        count= count+1;

%         plot(real((FDat(:,(fsiz(2)+1)-j,(fsiz(2)+1)-i,1,1,1,tf))),'LineWidth',1)


%         plot(real((FDat2(:,65-cr(i),tf))),'LineWidth',1)
        plot(real((FDat2(:,cr(i),tf))),'LineWidth',1)


        ylim([min(totDat(:)) max(totDat(:))])
%     end
end
set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
set(gca,'XTick',[], 'YTick', [])












