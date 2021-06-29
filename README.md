# k-DCE: k-means based Dynamic Contrast-Enhanced CT-scan segmentation

k-DCE performs k-means clustering for the segmentation of DCE-CT images. The program takes advantage of the differences in enhancing curves of the different tissues after the injection of contrast agent in order to automatically classify the image voxels. Clusters with specific enhancing curves corresponding to tissues are selected and the result is saved as the output segmentation image.  
The code is tuned to segment kidney cortex, medulla and tumor, but can also be used in other applications with minor adjustments about contrast values ranges. MRI data can also be used if pre-processed with histogram normalization.

## Dependencies 
	Matlab 2018b or higher


## Usage

Open the `main.m` script and set the parameters specific to your data. Compulsory variables to define are: 

	 img_folder   : the folder where you saved your Nifti images
	 t            : the acquisition times of your DCE data, in seconds

Run the `main.m` script in Matlab.