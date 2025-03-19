% Author:   Will Dong
% Date:     10/23/2024
% Version:  1.1
% Required packages: Image Processing Toolbox, Tools for NIfTI and ANALYZE
% image, Signal Processing Toolbox

% Purpose: To standardize the intensity values of all .nii files in myDir.

% Define the directory containing the NIfTI files
myDir = 'C:\Users\willd\Documents\Class_projects\Deep_learning\registration\unprocessed_cta\';

% Get a list of all *_id_original.nii files in the directory
niiFiles = dir(fullfile(myDir, '*_original.nii'));

% Loop through each file
for i = 1:numel(niiFiles)
    % Get the full file path
    CTIPath = fullfile(myDir, niiFiles(i).name);
    
    % Extract the ID from the filename
    [~, filename, ~] = fileparts(niiFiles(i).name);
    id = extractBefore(filename, '_original'); % Extract the ID from filename
    
    fprintf('Processing file: %s (ID: %s)\n', niiFiles(i).name, id);
    
    % Read the NIfTI image
    CTImage = niftiread(CTIPath);
    
    % Get voxel size from header
    headerInfo = niftiinfo(CTIPath);
    newHeaderInfo = headerInfo;
    
    % Add 1024 to all voxels in CTImage
    CTImage = CTImage + 1024;
    
    % Change all voxel values less than 0 to 10
    CTImage(CTImage < 0) = 10;
    
    % Save the corrected NIfTI file with the same naming convention
    outputFile = fullfile(myDir, ['corrected_' id '_original.nii']);
    niftiwrite(CTImage, outputFile, newHeaderInfo);
    
    fprintf('Saved corrected file: %s\n', outputFile);
end

fprintf('Processing complete for all files in %s\n', myDir);
