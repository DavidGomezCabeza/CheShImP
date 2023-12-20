from PIL import Image 
import numpy as np 
import matplotlib.pyplot as plt 
img = Image.open(r'D:/IBECPostDocDrive/2023_04_21_MatlabMRIDataProcessing/20230516_144623_Water_and_DMSO_eppendorf_MAR02002_DMSO_and_water_eppendorf_CSI_reconstruction_1_1/10/tmp_img/SpectraTimePoint_1_Phased.png') 
imgMRI = Image.open(r'D:/IBECPostDocDrive/2023_04_21_MatlabMRIDataProcessing/20230516_144623_Water_and_DMSO_eppendorf_MAR02002_DMSO_and_water_eppendorf_CSI_reconstruction_1_1/10/tmp_img/MRICut.png') 
img = img.convert("RGBA") 
datas = img.getdata() 
newData = [] 
for item in datas: 
    if item[0] == 255 and item[1] == 255 and item[2] == 255: 
        newData.append((255, 255, 255, 0)) 
    else: 
        newData.append(item) 
img.putdata(newData) 
plt.figure(figsize=(8,8)) 
plt.subplot(1, 1, 1) 
plt.imshow(imgMRI.resize([i*4 for i in np.shape(imgMRI)]), cmap='gray') 
plt.imshow(img.resize([i*4 for i in np.shape(imgMRI)])) 
plt.xticks([]) 
plt.yticks([]) 
plt.savefig(r'D:/IBECPostDocDrive/2023_04_21_MatlabMRIDataProcessing/20230516_144623_Water_and_DMSO_eppendorf_MAR02002_DMSO_and_water_eppendorf_CSI_reconstruction_1_1/10/tmp_img/MRICSIOverlay_1_Phased.png', bbox_inches='tight', pad_inches=0) 
