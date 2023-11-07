

function [namfil] = genSpectraAllColorGrid_Colormap_MaxRat(FDat, tf, foldpath, ppms)


    prompt = {'Min PPM Peak 1:','Max PPM Peak 1:', 'Min PPM Peak 2:','Max PPM Peak 2:'};
    dlgtitle = 'Ranges';
    dims = [1 35];
    definput = {'0','0','0','0'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    namfil = join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_Rat_',join(answer(:), '-'),'.png'], '');

    if ~isfile(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_Rat_',join(answer(:), '-'),'.fig'], ''))
        
        FDatTP = FDat(:,:,:,1,1,1,tf);
        
        FDatSiz = size(FDatTP);


        FDatTP = FDat(:,:,:,1,1,1,tf);

        poss = [length(ppms(ppms >= str2num(answer{2})))+1:length(ppms(ppms >= str2num(answer{1})))+1, length(ppms(ppms >= str2num(answer{4})))+1:length(ppms(ppms >= str2num(answer{3})))+1];
        allL = 1:length(ppms);

        for i = 1:length(allL)
            if ismember(i, poss)
                allL(i) = 0;
            end
        end

        nosMpre = FDatTP(allL(allL ~= 0 ),:,:);

        nosM = mean(nosMpre(:));
        nosS = std(nosMpre(:));

        pe1 = FDatTP(length(ppms(ppms >= str2num(answer{2})))+1:length(ppms(ppms >= str2num(answer{1})))+1,:,:);
        pe2 = FDatTP(length(ppms(ppms >= str2num(answer{4})))+1:length(ppms(ppms >= str2num(answer{3})))+1,:,:);


        rats = max(pe1)./max(pe2);

        FDatMaxs = max(pe1)./max(pe2);

        for i = 1:FDatSiz(2)
            for j = 1:FDatSiz(3)
                if (max(pe1(:,i,j))<= nosM+(nosS*4)) || (max(pe2(:,i,j))<= nosM+(nosS*4)) % or or and!?!?!?!?!
                    FDatMaxs(:,i,j) = 1;
                end
            end
        end

        FDatMaxsCmp = round(FDatMaxs(:)/min(FDatMaxs(:)));
        colmp = colormap(parula(length(unique(FDatMaxsCmp(:)))));
        colormap(gray)
%         colmp = colormap(autumn(length(unique(FDatMaxsCmp(:)))));
        sortDat = sort(unique(FDatMaxsCmp(:)));


        

        fsiz = size(FDat);
        totDat = real(FDat(:,:,:,1,1,1,:));
        totDatTP = real(FDat(:,:,:,1,1,1,tf));
    
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

%                     set(gca, 'color', [min([FDatMaxs1(1,j,i)/255, 1]), min([FDatMaxs2(1,j,i)/255, 1]), min([FDatMaxs3(1,j,i)/255, 1])]);
%                     set(gcf, 'color', [min([FDatMaxs1(1,j,i)/255, 1]), min([FDatMaxs2(1,j,i)/255, 1]), min([FDatMaxs3(1,j,i)/255, 1])]);    
                    mxd = round(max(FDatMaxs(:,j,i,1,1,1,tf))/min(FDatMaxs(:)));
                else % Rectangular grid
                    subplot(gensiz(3),gensiz(2),mm(count))
                    count= count+1;
                    plot(permute(real((FDat(:,i,j,1,1,1,tf))), [1,3,2,4,5,6,7,8]), 'black','LineWidth',1)

%                     set(gca, 'color', [min([FDatMaxs1(1,i,j)/255, 1]), min([FDatMaxs2(1,i,j)/255, 1]), min([FDatMaxs3(1,i,j)/255, 1])]);
%                     set(gcf, 'color', [min([FDatMaxs1(1,i,j)/255, 1]), min([FDatMaxs2(1,i,j)/255, 1]), min([FDatMaxs3(1,i,j)/255, 1])]); 
                    mxd = round(max(FDatMaxs(:,i,j,1,1,1,tf))/min(FDatMaxs(:)));
                end
                ylim([min(totDat(:))-(max(totDat(:))*0.20) max(totDat(:))])
                xlim([1-(gensiz(1)*0.10) gensiz(1)+(gensiz(1)*0.10)])


                
                set(gca, 'color', colmp(find(mxd == sortDat), :));
                set(gcf, 'color', colmp(find(mxd == sortDat), :)); 

                if FDatMaxs(:,i,j) == 1
                    set(gca, 'color', [0,0,0]);
                    set(gcf, 'color', [0,0,0]); 
                end

                
                set(gca,'XTick',[], 'YTick', [])
                set(gca, 'XColor','#e0dede', 'YColor','#e0dede')
            end
        end
        set(gca,'XTick',[], 'YTick', [])
        set(gcf, 'InvertHardcopy', 'off');
 

        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_Rat_',join(answer(:), '-'),'.png'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_Rat_',join(answer(:), '-'),'.fig'], ''))
        saveas(ff,join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'_Rat_',join(answer(:), '-'),'.svg'], ''))
        
        close(ff)
        
    else
%         openfig(join([foldpath,'\tmp_img\SpectraTimePoint_ColorGrid_',string(tf),'.fig'], ''),'visible');
    end

end


