# NTEC_Vessel_Analysis
Code for specific endothelial cells govern nanoparticle entry into solid tumours (Kingston et al., 2021). All code was evaluated with MATLAB 2019a with the DIP image toolbox with DIP image v2.7 (http://www.diplib.org or https://github.com/DIPlib/diplib).

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

# LICENSE

This software license is the 2-clause BSD license plus a third clause that prohibits redistribution and use for commercial purposes without further permission from the authors of this work (Kingston et al).

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Redistributions and use for commercial purposes are not permitted without the written permission of the authors of this work (Kingston et al). For purposes of this license, commercial purposes are the incorporation of the software into anything for which you will charge fees or other compensation or use of the software to perform a commercial service for a third party.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
