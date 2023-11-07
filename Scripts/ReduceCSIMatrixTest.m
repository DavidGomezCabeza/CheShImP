

%% Load Data

dat = load('C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\20230929_110110_DGC_HepG2_VolumeInitialTests_DGC_HepG2_VolumeInitialTests_1_1\3_CSI_20M\tmp_img\PhasedData.mat');
ppms = load('C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\20230929_110110_DGC_HepG2_VolumeInitialTests_DGC_HepG2_VolumeInitialTests_1_1\3_CSI_20M\tmp_img\ppms.mat').ppms;

tp=1;

magDat = dat.FDat3(:,:,:,1,1,1,tp);



%% Correct For Folding

% siz = size(magDat);
% 
% magDat = cat(2, magDat(:,2:end,:), magDat(:,1,:));


%% Zoom in Spectra

figure, 
plot(real(magDat(:,5,3)))

minin = 290;
maxin = 591;

magDat = magDat(minin:maxin, :, :);

ppms = ppms(minin:maxin);

%% Reduce CSI Matrix (sum Spectra)

xdi = 1;
ydi = 4;

sizZoom = size(magDat);

redDattmp = zeros(sizZoom(1), xdi, sizZoom(3));
redDat = zeros(sizZoom(1), xdi, ydi);

for i = 1:xdi
    
    in1 = siz(2)/xdi*i-siz(2)/xdi+1;
    in2 = siz(2)/xdi*i;

    redDattmp(:,i,:) = sum(magDat(:,in1:in2, :),2);

end

for i = 1:ydi
    
    in1 = siz(3)/ydi*i-siz(3)/ydi+1;
    in2 = siz(3)/ydi*i;

    redDat(:,:,i) = sum(redDattmp(:,:,in1:in2),3);

end


%% Normalise Spectra

for i = 1:xdi
    for j = 1:ydi

        normDat = (redDat(:,i,j)-min(redDat(:,i,j)))/(max(redDat(:,i,j))-min(redDat(:,i,j)));
        redDat(:,i,j) = normDat;
    end
end


%% Plot CSI

plotSpectraCSI(redDat, 1, '')


%% Extract Integral Ratios


lacMinInd = 77;
lacMaxInd = 100;

pyrMinInd = 191;
pyrMaxInd = 240;


% x = -2:0.001:2;
% 
% y = normpdf(x, 0, 0.2)
% 
% figure, plot(x, y)
% 
% A = trapz(x, y)


figure, plot(flip(ppms(lacMinInd:lacMaxInd)), flip(real(redDat(lacMinInd:lacMaxInd,i,j)))*1000)

figure, plot(flip(ppms(pyrMinInd:pyrMaxInd)), flip(real(redDat(pyrMinInd:pyrMaxInd,i,j)))*1000)


rats = zeros(xdi,ydi);
lacs = zeros(xdi,ydi);
pyrs = zeros(xdi,ydi);

mins = [];

for i = 1:xdi
    for j = 1:ydi
%         lacMax = trapz(flip(ppms(lacMinInd:lacMaxInd)), flip(real(redDat(lacMinInd:lacMaxInd,i,j))));
%         pyrMax = trapz(flip(ppms(pyrMinInd:pyrMaxInd)), flip(real(redDat(pyrMinInd:pyrMaxInd,i,j))));

        lacMax = trapz(flip(real(redDat(lacMinInd:lacMaxInd,i,j))));
        pyrMax = trapz(flip(real(redDat(pyrMinInd:pyrMaxInd,i,j))));

        rats(i,j) = lacMax/pyrMax;
        lacs(i,j) = lacMax;
        pyrs(i,j) = pyrMax;
        if lacMax < 0 
            mins = [mins, lacMax];
        end
        if pyrMax < 0 
            mins = [mins, pyrMax];
        end
    end
end

disp(rats')

lacs = lacs+(abs(min(mins)));
pyrs = pyrs+(abs(min(mins)));

disp((lacs./pyrs)')



%% Extract intensity ratios

lacMinInd = 85;
lacMaxInd = 109;

pyrMinInd = 200;
pyrMaxInd = 243;


rats = zeros(xdi,ydi);

for i = 1:xdi
    for j = 1:ydi
        lacMax = max(real(redDat(lacMinInd:lacMaxInd,i, j)));
        pyrMax = max(real(redDat(pyrMinInd:pyrMaxInd,i, j)));

        rats(i,j) = lacMax/pyrMax;
    end
end

disp(rats')





























