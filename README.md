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
  The small window will show the loaded proton image, and in green color the area or field of view (FoV) for the CSI selected. This is a non-interactive part of the platform. 

  ### 4.- 

  ### 5.- 

  ### 6.- 

  ### 7.- 

  ### 8.- 

  ### 9.- 

  ### 10.- 

  ### 11.- 

  ### 12.- 

  ### 13.- 

  ### 14.- 

  ### 15.- 

  ### 16.- 

  ### 17.- 

  ### 18.- 

  ### 19.- 

  ### 20.- 

  ### 21.- 

  ### 22.- 

  ### 23.- 

  ### 24.-

  ### 25.- 

  ### 26.- 

## References
  **1.- ds** sasa 
