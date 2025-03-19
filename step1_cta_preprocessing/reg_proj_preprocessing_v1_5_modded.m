% Author:   Will Dong (williamdong2029@u.northwestern.edu)
% Date:     10/23/2024
% Version:  1.5
% Required packages: Image Processing Toolbox, Tools for NIfTI and ANALYZE
% image, Signal Processing Toolbox

% Purpose:  To preprocess CT images in mlimage_0.nii format by resampling
% to a user-defined voxel size and cropping + zero-padding to a user-defined 
% field of view. 

% Inputs: image and ground truth mask nifti files. 
% Outputs: resampled and cropped image and ground truth mask nifti files,
% and excel sheet showing stats. 

% change directory to dataset location
%cd('C:\Users\willd\Documents\Markl_Lab\nu-cta-sample-data\preprocessing\')

% set directory to where the images and ground truth masks are present
myDir = 'C:\Users\willd\Documents\Class_projects\Deep_learning\registration\unprocessed_cta\intensity_standardized';

% initialize an empty cell array to store summary data for each image
outputData = {'Subject ID', 'Resampled Matrix Size x', 'y', 'z', ...
              'Cropped Matrix Size x', 'y', 'z', ...
              'Resampled Mask Volume (voxels)', 'Cropped Mask Volume (voxels)', ...
              'Cropped?', 'Inferior Crop?', 'Mask volume preserved?', ...
              'Lung peak position in upper half of image?', '% Volume saved', ...
              'Peaks above 60% lung threshold in lower half', 'MSE intensity difference in right vs left halves', ...
              'Peak mean intensity in upper half of image', 'posterior edge to crop distance (mm)', ...
              'anterior', 'right', 'left', 'superior', 'inferior', 'lung length in mm determined by 60% threshold', ...
              'Flagged for double checking'};

% iterate through files if it contains 'original.mlimage_0.nii' in its
% name. note that the 'original' represents only the CT image files, not
% the mask files. 

% this code is designed to work with files in the following naming
% convention: 

% idPart1_idPart2_original.mlimage_0.nii    image file
% idPart1_idPart2_mask.mlimage_0.nii        mask file

% e.g. 20100111_A8W10OK18_original.mlimage_0.nii for image file
% e.g. 20100111_A8W10OK18_mask.mlimage_0.nii for mask file

fileList = dir(fullfile(myDir, '*original.nii*'));
for k = 1:length(fileList)
baseFileName = fileList(k).name;
parts = split(baseFileName, '_');

% set id
id = strcat(parts{1}, '_', parts{2});

% declare voxel size to resample to (mm / voxel)
% note that in Matlab, x is right-left direction, y is anterior-posterior,
% and z is superior-inferior
desiredVoxelSize = [2.5, 2.5, 3];

% define crop distance for each direction, relative to the anchor point. 
% the anchor point is the middle slice for AP and RL axes and the upper 60% 
% threshold slice for SI axis (see diagram: last slice above threshold). 

cropSuperiorDistanceMM = 85;
cropInferiorDistanceMM = 435;
cropAnteriorDistanceMM = 100;
cropPosteriorDistanceMM = 155;
cropRightDistanceMM = 110;
cropLeftDistanceMM = 130;

% declare standardized (output) image matrix size (voxels). this can be
% determined by dividing the crop distance bounds (mm) by the voxel size
% (mm / voxel). then, a few voxels could be added to zero-pad the border,
% just in case rounding errors happen during resampling. 

% e.g.: 
% desiredVoxelSize = [2.5, 2.5, 3];
% cropSuperiorDistanceMM = 85;      z1
% cropInferiorDistanceMM = 435;     z2
% cropAnteriorDistanceMM = 100;     y1
% cropPosteriorDistanceMM = 155;    y2
% cropRightDistanceMM = 110;        x1
% cropLeftDistanceMM = 130;         x2

% [x, y, z] = [x1 + x2, y1 + y2, z1 + z2] ./ desiredVoxelSize + 2 desiredVoxelSize
% [101, 107, 180] = [110 + 130, 100 + 155, 85 + 435] ./ [2.5, 2.5, 3] + [5, 5, 6]

standardizedImageDims = [100, 104, 175];

%--------------------------------------------------------------------------

% open image and mask nifti files
CTImage = niftiread([id '_original.nii']);

% show volume rendering. This may be turned off to save time later on. 
%volshow(CTImage);
%volshow(CTMask);

% get voxel size from header
headerInfo = niftiinfo([id '_original.nii']);
CTIVoxelSize = headerInfo.PixelDimensions;


% resampling: 


% determine scaling factor
scalingFactors = CTIVoxelSize ./ desiredVoxelSize;
% calculate the new dimensions
newDims = round(size(CTImage) .* scalingFactors);
% resample the image and mask
resampledCTI = imresize3(CTImage, newDims);

% show volumes. Note that this matlab display does not account for new voxel
% size. To see results that are properly scaled, use 3D slicer. 
%volshow(resampledCTI);
%volshow(resampledCTM);

% declare new header information
newHeaderInfo = headerInfo;

newHeaderInfo.PixelDimensions = desiredVoxelSize;
newHeaderInfo.ImageSize = newDims;
newHeaderInfo.Datatype = 'uint16';

% save the resampled image as a new nifti
%niftiwrite(resampledCTI, [id '_resampled_CTI.nii'], newHeaderInfo);
%niftiwrite(resampledCTM, [id '_resampled_CTM.nii'], newHeaderInfoMask);

% define a body vs air threshold
greaterThanAirThrs = 300;

% generate binary image of the body
binaryImage = resampledCTI > greaterThanAirThrs;

%volshow(binaryImage)

% get the number of slices in image along the third dimension (axial)
numSlices = size(binaryImage, 3);

% get the last slice in the axial plane
lastAxialSlice = binaryImage(:, :, numSlices);


% display the last slice
%imshow(lastAxialSlice);
%title('Last Slice in the Axial Plane');

% fill in the trachea so the lung is closed. 
filledInTrachea = imfill(lastAxialSlice, "holes");

% display the last slice with trachea filled in
%imshow(filledInTrachea);
%title('Last Slice in the Axial Plane with trachea filled in');

% replace the last axial slice with the slice with filled in trachea 
binaryImage(:, :, numSlices) = filledInTrachea;

% fill in all the air pockets in the body
filledBinaryImage = imfill(binaryImage, "holes");

%volshow(filledBinaryImage)

% replace voxels in resampledCTI with 3000 where filledBinaryImage is true
% and the intensity in resampledCTI is less than the threshold
resampledAndFilledCTI = resampledCTI;
resampledAndFilledCTI(filledBinaryImage & (resampledCTI < greaterThanAirThrs)) = 3000;

%volshow(resampledAndFilledCTI);

% generate difference image, where all the air pockets are contrasted
differenceImg = resampledAndFilledCTI - resampledCTI;
% volshow(differenceImg);

% split into right and left difference image
numRLSlices = size(binaryImage, 2);
leftDifferenceImg = differenceImg(1:round(numRLSlices / 2), :, :);
rightDifferenceImg = differenceImg(round(numRLSlices / 2):end, :, :);

% calculate the mean intensity for each slice across dimensions 1 and 2
meanIntensity = mean(mean(differenceImg, 1), 2);
meanIntensityL = mean(mean(leftDifferenceImg, 1), 2);
meanIntensityR = mean(mean(rightDifferenceImg, 1), 2);

% convert to a vector (squeeze removes singleton dimensions)
meanIntensity = squeeze(meanIntensity);
meanIntensityL = squeeze(meanIntensityL);
meanIntensityR = squeeze(meanIntensityR);

% create a vector for the position (slice numbers)
numSlices = size(differenceImg, 3);
position = 1:numSlices;  % Assuming slices are numbered from 1 to numSlices

% calculate absolute difference 
intensityDifference = (meanIntensityL - meanIntensityR).^2;
MSE = mean(intensityDifference);

% Plot intensity vs position
%figure;
%plot(position, meanIntensity, 'LineWidth', 2);
%xlabel('Slice Position (Inferior-Superior Axis)');
%ylabel('Mean Intensity');
%title('Mean Intensity vs. Position');
%grid on;

% disregard the last slice, due to the filling in of the scanner table edge
meanIntensity(end) = 0;

% calculate the maximum value of the mean intensity
maxValue = max(meanIntensity);

% find slices where intensity is at least 60% of the maximum 
threshold = 0.6 * maxValue;
slicesAboveThreshold = find(meanIntensity >= threshold);

% Note to self: use findpeaks to see if there are multiple peaks. See the
% following link: https://www.mathworks.com/matlabcentral/answers/567648-how-to-detect-peaks-above-a-given-threshold-value
numOfPeaksAboveThreshold = length(findpeaks(meanIntensity(1:round(numSlices / 2)), "MinPeakHeight", threshold));

resampledMatrixSize = size(resampledCTI);        % Initial matrix size

% define first and last slice where the intensity is at least 60% of the
% max. this only works properly if there is one peak above the 60%
% threshold. 
firstSlice = slicesAboveThreshold(1);
lastSlice = slicesAboveThreshold(end);

%averageSlice = round((firstSlice + lastSlice) / 2);

sliceThicknessAP = desiredVoxelSize(2);
sliceThicknessRL = desiredVoxelSize(1);
sliceThicknessSI = desiredVoxelSize(3);

anteriorCropBound = round(resampledMatrixSize(2)/2 - cropAnteriorDistanceMM / sliceThicknessAP);
posteriorCropBound = round(resampledMatrixSize(2)/2 + cropPosteriorDistanceMM / sliceThicknessAP);
rightCropBound = round(resampledMatrixSize(1)/2 - cropRightDistanceMM / sliceThicknessRL);
leftCropBound = round(resampledMatrixSize(1)/2 + cropLeftDistanceMM / sliceThicknessRL);
superiorCropBound = round(lastSlice + cropSuperiorDistanceMM / sliceThicknessSI);
inferiorCropBound = round(lastSlice - cropInferiorDistanceMM / sliceThicknessSI);



% calculate the difference between the first and last slice
sliceRange = lastSlice - firstSlice;

% display the result
%disp(['First slice: ', num2str(firstSlice)]);
%disp(['Last slice: ', num2str(lastSlice)]);
%disp(['Difference between first and last slice: ', num2str(sliceRange)]);

% define lower bound position index to crop, based on 1.5 x lung size. 


%inferiorCrop = true;
if anteriorCropBound < 1
    anteriorCropBound = 1;
end
if posteriorCropBound > resampledMatrixSize(2)
    posteriorCropBound = resampledMatrixSize(2);
end
if rightCropBound < 1
    rightCropBound = 1;
end
if leftCropBound > resampledMatrixSize(1)
    leftCropBound = resampledMatrixSize(1);
end
if inferiorCropBound < 1
    inferiorCropBound = 1;
end
if superiorCropBound > resampledMatrixSize(3)
    superiorCropBound = resampledMatrixSize(3);
end


% crop image
croppedCTI = resampledCTI(rightCropBound:leftCropBound, anteriorCropBound:posteriorCropBound, inferiorCropBound:superiorCropBound);
%croppedCTM = imbinarize(croppedCTM);

%croppedCTI = resampledCTI(:, :, lowerCropPosition:end);
%croppedCTM = resampledCTM(:, :, lowerCropPosition:end);

%volshow(croppedCTI);

%volshow(resampledCTM);


% update dims

updatedDims = size(croppedCTI);

% create zerofilled images
zeroImage = zeros(standardizedImageDims, "uint16");
% calculate the starting indices for centering resampledCTI image onto the
% zerofilled image
startX = floor((standardizedImageDims(1) - updatedDims(1)) / 2) + 1;
startY = floor((standardizedImageDims(2) - updatedDims(2)) / 2) + 1;
startZ = floor((standardizedImageDims(3) - updatedDims(3)) / 2) + 1;

% paste the resampledCTI image into the center of the zero-filled image
zeroImage(startX:startX+updatedDims(1)-1, ...
          startY:startY+updatedDims(2)-1, ...
          startZ:startZ+updatedDims(3)-1) = croppedCTI;


% show volumes
%volshow(zeroImage);

% declare new header information for standardized image
newHeaderInfoFinal = newHeaderInfo;

newHeaderInfoFinal.PixelDimensions = desiredVoxelSize;
newHeaderInfoFinal.ImageSize = standardizedImageDims;

% save the standardized resampled image as a new nifti
niftiwrite(zeroImage, [id '_resampled_cropped_CTI_standardized.nii'], newHeaderInfoFinal);

% record values into table
resampledMatrixSize = size(resampledCTI);        % Initial matrix size
croppedMatrixSize = size(croppedCTI);     % Cropped matrix size





isImageCropped = true;

isMaskPreserved = true;

isLungPeakNormallyPositioned = true;
areThereLessThan10Peaks = true;
rLTotalDiffLessThan1000 = true;
%if totalDifference > 1000
%    rLTotalDiffLessThan1000 = false;
%end
lungLength = (lastSlice - firstSlice) * sliceThicknessAP;

% flagging to tell user to double-check
flagged = false;

outputData{end+1, 1} = id;                           % Subject ID
outputData{end, 2} = mat2str(resampledMatrixSize(1));     % Initial matrix size
outputData{end, 3} = mat2str(resampledMatrixSize(2));     % Initial matrix size
outputData{end, 4} = mat2str(resampledMatrixSize(3));     % Initial matrix size

outputData{end, 5} = mat2str(croppedMatrixSize(1));     % Cropped matrix size
outputData{end, 6} = mat2str(croppedMatrixSize(2));     % Cropped matrix size
outputData{end, 7} = mat2str(croppedMatrixSize(3));     % Cropped matrix size

outputData{end, 8} = 0;              % Initial mask volume
outputData{end, 9} = 0;              % Cropped mask volume
outputData{end, 10} = isImageCropped;                % T/F, is image cropped at all?
outputData{end, 11} = 'NA';                  % T/F, cropped in inferior direction?
outputData{end, 12} = isMaskPreserved;               % T/F, is mask volume preserved
outputData{end, 13} = isLungPeakNormallyPositioned;  % T/F, is lung positioned in upper half of image?
outputData{end, 14} = (prod(resampledMatrixSize) - prod(updatedDims)) / prod(resampledMatrixSize); % percent volume saved
outputData{end, 15} = numOfPeaksAboveThreshold; % Number of peaks above the 60% threshold located in lower half of image
outputData{end, 16} = MSE;       % MSE between right and left sides
outputData{end, 17} = maxValue; % peak mean intensity value in regions considered (upper 50% of image)
% end of for loop
outputData{end, 18} = 0; % Posterior edge distance
outputData{end, 19} = 0;  % Anterior edge distance
outputData{end, 20} = 0;     % Right edge distance
outputData{end, 21} = 0;      % Left edge distance
outputData{end, 22} = 0;  % Superior edge distance
outputData{end, 23} = 0;  % Inferior edge distance
outputData{end, 24} = lungLength;                  % Lung length in mm determined by 60% threshold
outputData{end, 25} = flagged;                  % flag to tell user to double check image (if there are peaks above threshold in the lower half of image or if the lung size is longer than 20 mm. 
end

writecell(outputData, 'CT_Preprocessing_Summary.xlsx');