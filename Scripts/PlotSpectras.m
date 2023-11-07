
figure
plot(ppms, real(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'b')
ylim([min(rd(:)) max(rd(:))])
xlim([min(ppms)+50 max(ppms)-40])
ax = gca;
ax.XDir = 'reverse';
xlabel('ppm')
legend('Real')


figure
plot(ppms,imag(handles.GUIDataAll.FDat2(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))), 'r')
ylim([min(id(:)) max(id(:))])
xlim([min(ppms)+50 max(ppms)-40])
ax = gca;
ax.XDir = 'reverse';
xlabel('ppm')
legend('Imag')


figure
plot(ppms+str2double(handles.CSC.String),handles.GUIDataAll.FDat(:,x2,y2,1,1,1,str2num(handles.TimePoints.String(handles.TimePoints.Value,:)))*str2double(handles.SFac.String), 'g')
ylim([min(fd(:)) max(fd(:))])
xlim([min(ppms+str2double(handles.CSC.String))+50 max(ppms+str2double(handles.CSC.String))-40])
ax = gca;
ax.XDir = 'reverse';
xlabel('ppm')
legend('Magn')