================================================================
PSBox --- A matlab toolbox for photometric stereo.
================================================================

Author: Ying Xiong.
Created: Jan 24, 2014.
Release: Jan 25, 2014 (v0.1).

================================================================
Quick start.
================================================================
>> addpath('NLLSBox');
>> addpath('Utils');
>> demoPSBox;

================================================================
Notation and convention.
================================================================

We use a right-hand coordinate system such that when an image is shown, x axis
points rightwards, y axis upwards and z axis outwards. The origin is therefore
at the bottom-left corner. The (x,y) pixel in an image 'I' is stored at
I(y,x). This convention only affects the I/O, where one needs to revert the y
axis, with following examples:
  # To display an image.
  imshow(I); axis xy;
  # To write an image to hard disk.
  imwrite(I(end:-1:1, :, :), 'filename.png');
  # To read an image from hard disk.
  I = im2double(imread('filename.png'));
  I = I(end:-1:1, :, :);

================================================================
Folder structure.
================================================================

TopDirectory
|-- Original               # Original images, will not be altered.
|
|-- OriginalRenamed        # Renamed original images.
|   |-- Image_NN.JPG       # Jpeg images.
|   |-- Image_NN.CR2       # Raw images, could be other suffix than 'CR2'.
|   |-- Image_NN.tiff      # Derendered raw images, could be other suffix
|   |                      # than 'tiff', e.g. 'png'. See 'rawOutSuffix'
|   |                      # variable in the code.
|   |-- ref.*              # Reference images.
|
+-- ManualData             # A folder containing manually extracted data.
|   |-- obj_bbox.txt       # The bounding box for the object.
|   |-- probes_bbox.txt    # The bounding box for light probes.
|   |-- circle{1,2}_pts.txt# The points on the circle of light probe.
|
+-- Objects                # Extracted object images.
|   |-- Image_NN.tiff      # Suffix can also be 'png' or others.
|
+-- LightProbe-{1,2}       # Extracted light probe images.
    |-- Image_NN.JPG
    |-- ref.JPG
    |-- circle_data.txt
    |-- light_directions.txt


================================================================
Features.
================================================================
* Fit the circle of chrome sphere from manual extracted points.
* Find lighting direction from given chrome sphere.
* Perform photometric stereo to recover albedo and normal map.
