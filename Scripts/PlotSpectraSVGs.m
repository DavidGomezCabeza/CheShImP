bw = imageObj2.Method.PVM_SpecSW(1);
bwc = imageObj2.Method.PVM_FrqWorkPpm(1);

   
ppms = flip(linspace(bwc-bw/2, bwc+bw/2, datsiz(1)));



figure
plot(ppms-5, FDat(:,x2,y2,1,1,1,1))
ax = gca;
xlim([155,205])

ax.XDir = 'reverse';
xlabel('ppm')








figure
plot(ppms, real(imageObj2.data(:,5,12,1,1,1,1)).^2 + imag(imageObj2.data(:,5,12,1,1,1,1)).^2)
ax = gca;
xlim([165,195])

ax.XDir = 'reverse';
xlabel('ppm')





ph2=0;
pivot = 211;
dsz = length(imageObj2.data(:,5,30,1,1,1,1));
figure
plot(ppms, real(imageObj2.data(:,5,30,1,1,1,1) .* exp(sqrt(-1).*((ccc(mada,1)+0.05+ph2.*((-pivot:-pivot+dsz-1)/dsz)'    )))))
ax = gca;
xlim([165,200])
ylim([min(real(imageObj2.data(:,5,30,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))-max(real(imageObj2.data(:,5,30,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))*0.05, ...
      max(real(imageObj2.data(:,5,30,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))*1])
ax.XDir = 'reverse';
xlabel('ppm')









figure
plot(ppms, real(imageObj2.data(:,x2,y2,1,1,1,1) * exp(sqrt(-1)*(-18))))
ax = gca;
xlim([165,220])

ax.XDir = 'reverse';
xlabel('ppm')




% 
% 
% figure
% plot(ppms, real(imageObj2.data(:,5,8,1,1,1,1) * exp(sqrt(-1)*(ccc(mada,1)+0.02))))
% ax = gca;
% xlim([165,190])
% ylim([min(real(imageObj2.data(:,5,8,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))-max(real(imageObj2.data(:,5,21,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))*0.005, ...
%       max(real(imageObj2.data(:,5,8,1,1,1,1) * exp(sqrt(-1)*ccc(mada,1))))*0.03])
% ax.XDir = 'reverse';
% xlabel('ppm')

