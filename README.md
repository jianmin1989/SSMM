# SSMM
A spatiotemporal shape-matching model for fusing MODIS and HLS data

	Zhang, X.; Wang, J.; Henebry, G.; Gao, F. 
	"Development and Evaluation of a New Algorithm for Detecting 30m Land Surface Phenology from VIIRS and HLS Time Series",
	ISPRS Journal of Photogrammetry and Remote Sensing.
	

	This code is developed by Dr. Jianmin Wang and Xiaoyang Zhang, South Dakota State University. 
	Copyright Dr. Xiaoyang Zhang and Jianmin Wang
 	Use of this code for commercial purposes is not permitted. 
	Contact Xiaoyang.Zhang@sdstate.edu or Jianmin.wang@sdstate.edu 
	for more information and updates of this code

	Version 1.0
	4 Jan 2020
	
----------------------------------------------------------------------------------------

The following files/directories should be present in directory HLSVII_fusion_release:
HLSVII_fusion_release - 
	HLS_VIIRS_shape_fusion.c
	HLS_VIIRS_shape_fusion.exe
	makefile
	input_example.txt
	input_description.txt
	example_data - (can be found in https://openprairie.sdstate.edu/global_land_surface_season_data/5/#attach_additional_files)
		HLS_EVI2_T18TYN_2016184.gz
		HLS_EVI2_T18TYN_2016187.gz
		...
		HLS_EVI2_T18TYN_2018181.gz
		VIIRS_EVI2.2016184.h12v04.BIP.gz
		VIIRS_EVI2.2016187.h12v04.BIP.gz
		...
		VIIRS_EVI2.2018181.h12v04.BIP.gz
			
			
HLSVII_fusion, Introductory Notes:
This application is developed to fuse HLS and VIIRS (or MODIS) vegetation index (VI) time series.  


USAGE AND INPUTS
Usually, you can directly run "HLS_VIIRS_shape_fusion.exe input.txt" to fuse HLS and VIIRS VI time series, after providing the inputs in "input.txt". The users can also modify the c code "HLS_VIIRS_shape_fusion.c" and compile it with the makefile to update the "HLS_VIIRS_shape_fusion.exe". 
	 
On running, it takes the file 'input.txt' to specify the 13 input parameters. "input_example.txt" is provided as an example of input.txt. The explaination of required parameters can be found in "input_description.txt." 
			
All the input (HLS and VIIRS) and output images should be stored as gz files. The length of VI time series (number of images) is the same for HLS and VIIRS. The input HLS and VIIRS images must be stored in separate files for different dates. The dates should be provided as the 13th parameter. HLS and VIIRS VI file names should include the date in the format the same to the 13th parameter. In the file name templets of HLS and VIIRS, date should be indicated as "DDDDDDD". 
VI values of HLS and VIIRS images should be scaled by 10,000 and stored as 16-bits integer (from -32768 to 32767). 

For HLS images, the map coordinates of upperleft corner and UTM zone must be provided to obtain the projection and mapping information. For VIIRS images, the projection and mapping information will be obtained directly from the tile name. 
One HLS tile and all the VIIRS tiles covering that HLS tile should be provided. One HLS can be covered by at most 4 VIIRS tiles. The users should provide the vertical and horizontal number of available VIIRS tiles in the input parameters. Only the pixels covered in HLS tile covered by the available VIIRS tiles will be fused. In the VIIRS file name templet, tile name should be indicated as "hHHvVV".

 


