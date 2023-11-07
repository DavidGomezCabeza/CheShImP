%% Function to process the FID data with basic automatic phase correction

function [imageObj3, sizz, epc, ephci, siss, phv, pap, phma, dsz, ph1, pivotppm, pivot, handles] = processPhaseDat(kdataObj2, imageObj2, handles)


    imageObj3=kdataObj2.reco('quadrature');
    sizz = size(imageObj3.data);
    epc = [];
    for i = 1:sizz(2)
        for j=1:sizz(3)
            mm=imageObj3.data(:,i,j,1,1,1,1);
            df = abs(diff(abs(real(mm))));
            th = mean(df);
            ms = find(df'>(th-4*std(th))');
            epc = [epc, ms(1)];
        end
    end
    
    ephci = floor(mean(epc));


    imageObj3=imageObj3.reco('phase_rotate');
    imageObj3=imageObj3.reco('zero_filling');

    % THIS IS FOR PHASE CORRECTION!!! You have to remove phase_corr_pi!
    siss = size(imageObj3.data);
    for i = 1:siss(2)
        for j = 1:siss(3)
            if length(siss) >3
                for z = 1:siss(7)
                    imageObj3.data(:,i,j,1,1,1,z) = [imageObj3.data(ephci+1:end,i,j,1,1,1,z); imageObj3.data(1:ephci,i,j,1,1,1,z)];
%                         kdataObj2.data(:,i,j,1,1,1,z) = [kdataObj2.data(ephci+1:end,i,j,1,1,1,z); kdataObj2.data(1:ephci,i,j,1,1,1,z)];
                end
            else
                imageObj3.data(:,i,j) = [imageObj3.data(ephci+1:end,i,j); imageObj3.data(1:ephci,i,j)];
%                     kdataObj2.data(:,i,j) = [kdataObj2.data(ephci+1:end,i,j); kdataObj2.data(1:ephci,i,j)];
            end
        end
    end
    


    imageObj3.data = fft(imageObj3.data);
    imageObj3=imageObj3.reco('cutoff');
    imageObj3=imageObj3.reco('scale_phase_channels');
    imageObj3=imageObj3.reco('transposition');
    imageObj3.data = flip(imageObj3.data,1);


    phv = -8:0.1:8;
    siss = size(imageObj3.data);

    if length(siss) > 3
        pap = siss(7);
    else
        pap = 1;
    end

    phma = zeros(siss(2), siss(3), pap);


    dsz = length(imageObj2.data(:,1,1,1,1,1,1));
    ph1 = zeros(siss(2), siss(3), pap);
    pivotppm = zeros(siss(2), siss(3), pap) + str2num(handles.PivotEdit.String);
    pivot = zeros(siss(2), siss(3), pap);


    for z = 1:pap
        for i = 1:siss(2)
            for j = 1:siss(3)

                    ccc = zeros(length(phv), 2);
                    ccc(:,1) = phv;
                    for k =1:length(phv)
%                         phdat = real(imageObj3.data(:,i,j,1,1,1,z) .* exp( ...
%                                 sqrt(-1) .* (phma(i,j,z) + ph1(i,j,z) .* ( ...
%                                 (-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)')));
                        phdat = real(imageObj3.data(:,i,j,1,1,1,z) .* exp(sqrt(-1)*phv(k) +ph1(i,j,z).*((-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)'   ));
                        
                        ccc(k,2) = sum(phdat(phdat<0));
                
                    end
                    mada= find(ccc(:,2) == max(ccc(:,2)));

                    phv2 = phv(mada)-0.1:0.001:phv(mada)+0.1;

                    ccc2 = zeros(length(phv2), 2);
                    ccc2(:,1) = phv2;
                    for k =1:length(phv2)
%                         phdat = real(imageObj3.data(:,i,j,1,1,1,z) .* exp( ...
%                                 sqrt(-1) .* (phma(i,j,z) + ph1(i,j,z) .* ( ...
%                                 (-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)')));
                        phdat = real(imageObj3.data(:,i,j,1,1,1,z) .* exp(sqrt(-1)*phv2(k) +ph1(i,j,z).*((-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)' ));
                        
                        ccc2(k,2) = sum(phdat(phdat<0));
                
                    end
                    mada2= find(ccc2(:,2) == max(ccc2(:,2)));

                    phma(i,j,z) = ccc2(mada2(1),1);
                    imageObj3.data(:,i,j,1,1,1,z) = imageObj3.data(:,i,j,1,1,1,z) .*exp(sqrt(-1)*ccc2(mada2(1),1) +ph1(i,j,z).*((-pivot(i,j,z):-pivot(i,j,z)+dsz-1)/dsz)' );
            end
        end
    end


    % Base Line Correction
%     blp = abs(real(max(kdataObj2.data(:))))/2;
%     blp = 50000;
%     for z = 1:pap
%     %         disp(z)
%         for i = 1:siss(2)
%             for j = 1:siss(3)
%                 kdataObj2.data(:,i,j,1,1,1,z) = kdataObj2.data(:,i,j,1,1,1,z)/blp;
%                 imageObj.data(:,i,j,1,1,1,z) = imageObj.data(:,i,j,1,1,1,z)/blp;
%                 imageObj2.data(:,i,j,1,1,1,z) = imageObj2.data(:,i,j,1,1,1,z)/blp;
%                 imageObj3.data(:,i,j,1,1,1,z) = imageObj3.data(:,i,j,1,1,1,z)/blp;
%             end
%         end
%     end





end