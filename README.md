# NTEC_Vessel_Analysis
Specific endothelial cells govern nanoparticle entry into solid tumours (Kingston et al., 2021)

# Convert .czi files to multipage .tiff files for each imaging channel
This is done using the convert_czi_to_ometiff_revised_func.m fuction. 

# Pre-processing the blood vessel channel
A local and global intensity normalization is done throughout the imaging volume to prepare the blood vessel channel for segmentation using Ilastik. This is done using the pre_process_vess_func.m fuction. 

# Blood vessel segmentation 
The blood vessel segmenation was done using the pixel classification tool in Ilastik (https://www.ilastik.org/ or https://github.com/ilastik/). 

# Blood vessel analysis
The segmented blood vessels and nanoparticle channels are used to define N-TEC locations along the blood vessels and quantify the blood vessel morphology and local nanoparticle intensity using the vessel_ntec_analysis_func.m function. 

# Simulation of nanoparticle access from N-TECs/tumour blood vessels 
Access to the tumour space is simulated using 3D images of tumour blood vessels using the distributionfromntecs_func.m fuction. Three simulations are done: 1) nanoparticle access from actual N-TEC locations, 2) nanoparticle access from random N-TEC locations, and 3) nanoparticle access from all vessel surfaces. 
