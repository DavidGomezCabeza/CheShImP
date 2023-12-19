# CheShImP

![alt text](https://github.com/DavidGomezCabeza/CheShImP/blob/main/Logo.png?raw=true)

CheShImP is a program to process Chemical Shift Images (CSI) from Magnetic Resonance Imaging (MRI) systems. The platform automatically processes the FID files (for now, only tested with Bruker systems using ParaVision 360) and reconstructs the spectroscopic image for the users. The program can also overlay these CSIs over an MRI proton image (DICOM format), and the platform will automatically save all images and overlays. By default, the platform works with magnitudes by default, but manual/automatic phase correction is also available. 

The system has been tested for Bruker systems using Paravision 360, reading the fid_proc file. Nonetheless, we are working on expanding it to different versions and systems. 

## Installation

For non matlab users, the program (and its GUI) can be installed using the [**CheShImPInstaller_web.exe**](https://github.com/DavidGomezCabeza/CheShImP/tree/main/CheShImP/for_redistribution) executable file intended for redistribution. 

Matlab users can access to the source code located in the directory [**Scripts**](https://github.com/DavidGomezCabeza/CheShImP/tree/main/Scripts). In there, executing the main script **ImageMRIDat_DVD.m** starts the CheShImP GUI. 

**IMPORTANT**

To generate the image overlay between CSI and a proton image, the platform uses python. For this reason, it is important to have python installed in your PATH System Variables. This Python Installation needs to have installed PIL, numpy and matplotlib. 

## Package Functionalities

![alt text](https://github.com/DavidGomezCabeza/CheShImP/blob/main/Platform.png?raw=true)

  ### 1.- CSI Path Directory. 
  Select the MAIN directory path for your CSI data (usually named with a number in ParaVision 360 for Bruker systems). The button will open a window search and once inside the directory, just select a file from it (does not matter which one) to load you data and start the processing. After this, you will be able to use the rest of features of the platform. 

  ### 2.- Proton Image Path Directory
  Select the MAIN directory path for your proton image to overlay under your CSI. The button will open a window search and once inside the directory, just select a file from it (does not matter which one) or select a specific DICOM file (recomended option). The platform will read the Bruker generated meta-data files, so it will generate the overlay with the correct FOV from the sample (hence, the CSI field of view needs to be define inside the area of the proton image used). 

  ### 3.- CSI and Proton Image FoV Overlay
  The small window will show the loaded proton image, and in green colour the area or field of view (FoV) for the CSI selected. This is a non-interactive part of the platform. 

  ### 4.- Time-Point Selection
  If you selected more than one repetition for the spectroscopic image in your MRI system, the platform will allow you to select the desired one to display. 
  
  ### 5.- Voxel Selected Position 
  Whenever you select a specific voxel in the grid displayed in panel 13, this part of the platform will display the X and Y coordinates for it to help the user navigate the platform. 

  ### 6.- SNR Computation and Display of the Best One
  This section of the platform displays the highest signal-to-noise ratio (SNR) of the image (i.e. the one for the voxel with the highest one). Nonetheless, the platform will save a Matlab file named SNR.mat inside the tmp_img directory with the SNR for each individual voxel. To automatically compute the SNR, the platform takes 10% of each edge from an individual spectrum to extract the standard deviation of the noise, and uses the intensity of the highest peak as the signal. If the user wants to reduce this edge threshold, it can be modified in the scrip compSNR.m.

  ### 7.- Overlay of CSI and Proton Image


  ### 8.- Scaling Factor (zoom)


  ### 9.- Line Broadening


  ### 10.- Chemical Shift Correction


  ### 11.- Switch Between Magnitude and Phase Mode


  ### 12.- Phase (0 and 1) and base-line correction


  ### 13.- CSI Display and Voxel Selection


  ### 14.- CSI Colour Grid Displays (non-interactive)


  ### 15.- Real Component Plot for Selected Voxel (fft as default)


  ### 16.- Imaginary Component Plot for Selected Voxel (fft as default)


  ### 17.- Magnitude or Phased Spectra Plot for Selected Voxel


  ### 18.- Restart GUI for New Experiment


  ### 19.- Delete Temporary Directory with Processed Data


  ### 20.- Plot FID of Selected Voxel (Displayed in Windows for Points 15 and 16)


  ### 21.- Use Colour Grid Instead of Porton Image in Overlay for Window of Point 13


  ### 22.- Generate Colour Grid Using Maximum Peak Intensity in Each Voxel


  ### 23.- Generate RGB Colour Grid Using Maximum of a Spectral Region for Each Primary Colour


  ### 24.- Level of Saturation Applied in the RGB Colour Grid (See section 23)


  ### 25.- Peak Intensity Value Selected in Window for Point 17


  ### 26.- Chemical Shift in Parts Per Million (ppm) Value Selected in Window for Point 17


## References
  **1.- ds** sasa 
