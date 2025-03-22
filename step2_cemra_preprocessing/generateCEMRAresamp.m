

% open image and mask nifti files
%CTImage = niftiread("C:\Users\willd\Documents\Markl_Lab\nu-cta-sample-data\preprocessing\corrected_613_resampled_cropped_CTI_standardized.nii");
CEMRAImage = niftiread("C:\Users\willd\Downloads\COR_CEMRA_TT_1_0s.nii.gz");

% show volume rendering. This may be turned off to save time later on. 
%volshow(CTImage);
%volshow(CTMask);

% get voxel size from header
headerInfo = niftiinfo("C:\Users\willd\Downloads\COR_CEMRA_TT_1_0s.nii.gz");
CEMRAVoxelSize = headerInfo.PixelDimensions;

% resampling: 
desiredVoxelSize = [2.5, 2.5, 3];


% determine scaling factor
scalingFactors = CEMRAVoxelSize ./ desiredVoxelSize;
% calculate the new dimensions
newDims = round(size(CEMRAImage) .* scalingFactors);
% resample the image and mask
resampledCEMRA = imresize3(CEMRAImage, newDims);

headerInfo.ImageSize = size(resampledCEMRA);
headerInfo.PixelDimensions = desiredVoxelSize;

niftiwrite(resampledCEMRA, 'TESTresampled_CEMRA.nii', headerInfo);
