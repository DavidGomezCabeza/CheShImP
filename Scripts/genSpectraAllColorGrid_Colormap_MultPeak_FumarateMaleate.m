

function [namfil] = genSpectraAllColorGrid_Colormap_MultPeak_FumarateMaleate(FDat, tf, foldpath, ppms)


    prompt = {'R- Min PPM Peak 1:','R- Max PPM Peak 1:', 'G- Min PPM Peak 2:','G- Max PPM Peak 2:', 'B- Min PPM Peak 3:','B- Max PPM Peak 3:'};
    dlgtitle = 'Ranges';
    dims = [1 35];
    definput = {'0','0','0','0','0','0'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    namfil = join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_RGB_',join(answer(:), '-'),'.png'], '');

    if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_RGB_',join(answer(:), '-'),'.fig'], ''))
        
        sat = 0.25;

        FDatTP = FDat(:,:,:,1,1,1,tf);
        
        FDatSiz = size(FDatTP);

        if str2num(answer{1}) == 0 && str2num(answer{2}) == 0
            FDatMaxs1 = zeros(1,FDatSiz(2), FDatSiz(3));
        else
            xx = max(FDatTP(length(ppms(ppms >= str2num(answer{2})))+1:length(ppms(ppms >= str2num(answer{1})))+1,:,:));

            FDatMaxs1 = round(255*(xx)/(max(xx(:))*sat));
            
        end

        if str2num(answer{3}) == 0 && str2num(answer{4}) == 0
            FDatMaxs2 = zeros(1,FDatSiz(2), FDatSiz(3));
        else
            yy = max(FDatTP(length(ppms(ppms >= str2num(answer{4})))+1:length(ppms(ppms >= str2num(answer{3})))+1,:,:));
            zz = max(FDatTP(length(ppms(ppms >= str2num(answer{6})))+1:length(ppms(ppms >= str2num(answer{5})))+1,:,:));
            
            zit = size(yy);
            for j = 1:zit(2)
                for k = 1:zit(3)
                    if yy(1,j,k) < 0
                        yy(1,j,k) = 0;
                    end
                    if zz(1,j,k) < 0
                        zz(1,j,k) = 0;
                    end
                end
            end

            dd= yy+zz;
            FDatMaxs3 = round(255*(dd)/(max(dd(:))*sat));
        end

        FDatMaxs2 = zeros(1,FDatSiz(2), FDatSiz(3));
        

        colmp = [1:255; 1:255; 1:255]'/255;
        colormap(gray)






        fsiz = size(FDat);
        totDat = real(FDat(:,:,:,1,1,1,:));
        totDatTP = real(FDat(:,:,:,1,1,1,tf));
    
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
                if gensiz(2)==gensiz(3) % squared grid
                    subplot(gensiz(2),gensiz(3),count)
                    count= count+1;
                    plot((real((FDat(:,j,i,1,1,1,tf)))), 'black','LineWidth',1)

                    set(gca, 'color', [min([FDatMaxs1(1,j,i)/255, 1]), min([FDatMaxs2(1,j,i)/255, 1]), min([FDatMaxs3(1,j,i)/255, 1])]);
                    set(gcf, 'color', [min([FDatMaxs1(1,j,i)/255, 1]), min([FDatMaxs2(1,j,i)/255, 1]), min([FDatMaxs3(1,j,i)/255, 1])]);    

                else % Rectangular grid
                    subplot(gensiz(3),gensiz(2),mm(count))
                    count= count+1;
                    plot(permute(real((FDat(:,i,j,1,1,1,tf))), [1,3,2,4,5,6,7,8]), 'black','LineWidth',1)

                    set(gca, 'color', [min([FDatMaxs1(1,i,j)/255, 1]), min([FDatMaxs2(1,i,j)/255, 1]), min([FDatMaxs3(1,i,j)/255, 1])]);
                    set(gcf, 'color', [min([FDatMaxs1(1,i,j)/255, 1]), min([FDatMaxs2(1,i,j)/255, 1]), min([FDatMaxs3(1,i,j)/255, 1])]); 

                end
                ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
                xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])


                

                
                set(gca,'XTick',[], 'YTick', [])
                set(gca, 'XColor','#e0dede', 'YColor','#e0dede')
            end
        end
        set(gca,'XTick',[], 'YTick', [])
        set(gcf, 'InvertHardcopy', 'off');
 

        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGridMalate_',string(tf),'_RGB_',join(answer(:), '-'),'.png'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGridMalate_',string(tf),'_RGB_',join(answer(:), '-'),'.fig'], ''))
        
        close(ff)
        
    else
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');
    end

end


