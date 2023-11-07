
x1 = data(:,:,:,:,channel,NI,NR);
x2 = reco_phase_rotate( data(:,:,:,:,channel,NI,NR) , Reco, map(NI, NR) );

figure
hold on
plot((x1(:,3,3)))
plot((x2(:,3,3)))


figure
hold on
plot(real(x1(:,3,3)))
plot(real(x2(:,3,3)))





