







i=4;
j=4;
z=1;

epc = [];

for i = 1:8
    for j=1:8

    mm=imageObj3.data(:,i,j,1,1,1,1);

% figure, plot(real(mm))
% figure, plot(abs(real(mm)))
% figure, plot(diff(abs(real(mm))))



df = abs(diff(abs(real(mm))));

th = mean(df);



% figure, plot(df(df>th))




% df'>(th-3*std(th))'


ms = find(df'>(th-4*std(th))');

% disp(ms(1))

epc = [epc, ms(1)];

    end
end

floor(mean(epc))




