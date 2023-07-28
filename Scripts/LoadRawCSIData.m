
fileID = fopen(join(['C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\UCSFtest_processed','\fid'],''),'r');
fidFile=fread(fileID, 'int32');
fclose(fileID);



fileID2 = fopen(join(['C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\UCSFtest_processed\dyn_fid_01','\sfid'],''),'r');
fidFile2=fread(fileID2, 'int32');
fclose(fileID2);


fileID3 = fopen(join(['C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\UCSFtest_processed\dyn_fid_01','\ser'],''),'r');
fidFile3=fread(fileID3, 'int32');
fclose(fileID3);



fileID = fopen(join(['C:\IBECPostDocDrive\2023_04_21_MatlabMRIDataProcessing\VoxelsPossitionTests\20230426_170750_CSIVoxelTests_2_CSIVoxelTests_2_1_1\5\','\fid_proc.64'],''),'r');
fidFile=fread(fileID, 'int32');
fclose(fileID);

fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));


FDat = reshape(fidFile2, 128*2,8,8)

figure, 
plot(fidFile2(129:128*2))















fileID = fopen(join(['D:\20230510_170210_jevo_calibration_VOL13C_calibration3_DMSO_BrukerPhantom_vol_1_1\7','\rawdata.job0'],''),'r');
fidFile=fread(fileID, 'int32');
fclose(fileID);

fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));

fidFile3 = reshape(fidFile2, 1, 2048, 16*16);




fileID = fopen(join(['D:\20230510_170210_jevo_calibration_VOL13C_calibration3_DMSO_BrukerPhantom_vol_1_1\7','\fid_proc.64'],''),'r');
fidFile=fread(fileID, 'int64');
fclose(fileID);

fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));

fidFile3 = reshape(fidFile2, 1, 2048, 16*16);








% nPoints x nRx x nSlices x nAv x nPh1 x nPh2 x nRep

fileID = fopen(join(['D:\20230510_170210_jevo_calibration_VOL13C_calibration3_DMSO_BrukerPhantom_vol_1_1\13','\rawdata.job0'],''),'r');
fidFile=fread(fileID, 'int32');
fclose(fileID);

fidFile(10000)

fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));

fidFile3 = reshape(fidFile2, 1, 2048, 16*16);








fileID = fopen(join(['D:\20230510_170210_jevo_calibration_VOL13C_calibration3_DMSO_BrukerPhantom_vol_1_1\7','\2dseq'],''),'r');
fidFile=fread(fileID, 'int8');
fclose(fileID);

fidFile2=complex(fidFile(1:2:end), fidFile(2:2:end));

fidFile3 = reshape(fidFile2, 1, 2048, 16*16);















