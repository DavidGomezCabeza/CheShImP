
%% This is for fumarate
ms = size(FDat);

%% SNR calculation
FDat2 = reshape(FDat, ms(1), ms(2)*ms(3));

SNRs = zeros(1,ms(2)*ms(3));
for i = 1:ms(2)*ms(3)

    noi = [FDat2(1:150,i); FDat2(350:end,i)]; % This is for fumarate
    snr = max(FDat2(:,i))/std(noi);
    SNRs(i) = snr;

end


%% Peak aplitude calculation


for i = 1:ms(2)*ms(3)

    noi = [FDat2(1:150,i); FDat2(350:end,i)];

    mean(noi)

    snr = max(FDat2(:,i))/std(noi);
    SNRs(i) = snr;

end


figure, scatter(ppms, FDat2(:,201))
hold on
plot([ppms(1) ppms(end)],[(max(FDat2(:,i))-mean(noi))*0.5 (max(FDat2(:,i))-mean(noi))*0.5])
plot([ppms(1) ppms(end)],[(max(FDat2(:,i))-mean(noi))*0.1 (max(FDat2(:,i))-mean(noi))*0.1])


diff(ppms(1:4))








