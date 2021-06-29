% Main script to be launched in order to segment the 
% Dynamic contrast-Enhanced CT-scan images (DCE-CT) by voxels clustering
% Images must be registered and with same resolution
% (ANTs does this automatically)

clear variables
close all

%% Define your data variables. 

% Define the folder with Niftifiles to be segmented
% Example:
img_folder = 'DATA/';

% Set acquisition times of different phases (in seconds)
% Example:
t = [0  237  291  878]; 


%% --- Parameters for image pre-processing ---

% Thresholds of values to be ignored
% Default values tuned on kidney CT-scan:
t_high = 300;
t_low = -50;

% Set to 0 to load images already computed, to speed up the code
dif_avgnbr = 1; % Computation of the difference between a voxel and the avg of its neighbours
Im_all_preproc = 1; % Preprocess images to remove non interesting voxels

% Set to 0 to skip plotting of preprocessed images
plotting = 1;

%% Load Nifti images 
load_registered_nii


%% Apply clustering 
% Set number of clusters (otherwise it will be automatically estimated)
n_clusters = 6; 
clustering3D


%% save Nifti image with clusters labels
niftiwrite(single(pixel_labels), [folder_out 'segmentation.nii'], info)