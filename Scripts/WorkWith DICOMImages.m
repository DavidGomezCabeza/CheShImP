
%%
dicominfo("E:\IBEC\NAFLD_PostDoc\2023_04_21_MatlabMRIDataProcessing\MRIImages\MRIm3.dcm")
[X,map] = dicomread("E:\IBEC\NAFLD_PostDoc\2023_04_21_MatlabMRIDataProcessing\MRIImages\MRIm3.dcm");

montage(X,map,"Size",[1 1]);



%%
info = dicominfo("E:\IBEC\NAFLD_PostDoc\2023_04_21_MatlabMRIDataProcessing\MRIImages\MRIm3.dcm");
Y = dicomread(info);
dim = size(Y);
figure, 
imshow(Y,[])
hold on
scatter(dim(1)/2, dim(2)/2)











FOV = header.PVM_Fov;
nPhaseEncodes = header.PVM_EncMatrix;
xAxis = linspace(-FOV(1)/2,FOV(1)/2,nPhaseEncodes(1));
yAxis = linspace(-FOV(2)/2,FOV(2)/2,nPhaseEncodes(2));

xAxisPix = dim(1)*xAxis/FOV(1)+128
yAxisPix = dim(2)*yAxis/FOV(2)+128
yAxisPix=yAxisPix*2

figure
grida = 0:8;
grid1 = [xAxisPix;xAxisPix];
grid2 = repmat([xAxisPix(1);xAxisPix(end)],1,length(xAxisPix));
% This seems wuite promissing
% figure

figure, 
% imshow(Y,[])
hold on
plot([xAxisPix;xAxisPix],repmat([yAxisPix(1);yAxisPix(end)],1,length(yAxisPix)),'g','LineWidth',2)
plot(repmat([xAxisPix(1);xAxisPix(end)],1,length(xAxisPix)),[yAxisPix;yAxisPix],'g','LineWidth',2)
% xlim([0 256])
% ylim([0 256])

ff=rectangle('Position',[73.1429 73.1429 73.1429/2 73.1429], 'FaceColor','g', 'EdgeColor', 'none') %  [x y w h]

x = [7 8 8 7];
y = [7 7 8 8];
figure, fs=patch(x,y,'g', 'EdgeColor', 'none','FaceAlpha',.2)


figure
grida = 0:8;
grid1 = [grida;grida];
grid2 = repmat([0;8],1,length(grida));
% This seems wuite promissing
% figure
hold on
plot(grid1,grid2,'g','LineWidth',2)
plot(grid2,grid1,'g','LineWidth',2)




























