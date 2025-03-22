# CEMRA-to-CTA
 
Instructions: 

Step 1: Preprocess the CTA images

Run the reg_proj_preprocessing_v1_5_modded.m script in Matlab to automatically resample and crop your CTA images to a voxel size of 2.5 mm x 2.5 mm x 3 mm and a matrix size of 100 x 104 x 175. 


Step 2: Preprocess the CEMRA images

Run the generateCEMRAresamp.m script in Matlab to resample your CEMRA images to a voxel size of 2.5 mm x 2.5 mm x 3 mm. 


Step 2.5: Registration and intensity standardization

Open 3D slicer and use General Registration (ANTs) module to register the CEMRA image to the corresponding CTA image. The CTA image is the fixed image and CEMRA image is the moving image. Use the deformable registration (SyN) and set inital transform as image center of mass. Newer version of 3D slicer may be bugged in running ANTs, but version 5.6.2 is known to work. 

For intensity standardization, set the intensity range of all images to the [-1, 1] range. 


Step 3: Post registration processing

Offical method: Use a Python script to crop out 15 sagittal slices on all CTA and CEMRA images. This should make the matrix size 70 x 104 x 175. 

Experimental method: Use prototypes\post_registration_processing_intersection_mrthod.ipynb to zero fill any regions of CTA and CEMRA that doesn't intersect. Remove all blank slices. 

Step 4: Train and evaluate model

Run monai_model_training_and_eval. Consider transferring the code to be run in a high performance computing cluster environment, as it is quite computationally expensive. A GPU would be necessary. 

Step 5: Post model processing

Run post_model_processing. This will recover the voxel size lost during training. 
