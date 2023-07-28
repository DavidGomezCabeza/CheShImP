
bruker_version % To check the version of the package

edit BrukerExample % To open the main example script. All quite well explained in there. See manual for more details. 

% Reconstructing imported k-space data section does not work with the
% example? I think that the code has been touched, so it would be good to
% get in touch with the lab and get the original code to compare. The issue is that the reco object is always empty.

fileattrib E:\IBEC\NAFLD_PostDoc\2023_04_21_MatlabMRIDataProcessing\pvmatlab\export\Kiwi_pvmatlab\2\pdata\1\2dseq % This is to check the file permissions


fileattrib "pvmatlab\export\Kiwi_pvmatlab\2\pdata\1\2dseq"



kdataObj = kdataObj.setDataPath('reco','TestData/2/pdata/1/reco');
kdataObj = kdataObj.readReco;
imageObj = kdataObj.reco('all','image');



imageObj = imageObj.genExportVisu('genmode','Visu', ...
['TemplatePath', 'E:\IBEC\NAFLD_PostDoc\2023_04_21_MatlabMRIDataProcessing\UCSFtest_nonprocessed\visu_pars']);


figure, plot(dataset(:,1,1,1,1,1,1))

figure, plot(real(dataset(:,1,1,1,1,1,1))), hold on
plot(imag(dataset(:,1,1,1,1,1,1)))











