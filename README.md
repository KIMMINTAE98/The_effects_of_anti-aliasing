# The effects of anti-aliasing
This program shows effects of anti-aliasing about YUV 4:2:0 8bit format images.

Input  : original image

Output : 1/2 downsampled image with/without anti-aliasing, reconstructed image with/without anti-aliasing

----------------

## Usage
set input image path
``` C++
string input_file = "./Input/Cactus_1920x1080_yuv420_8bit_frame200.yuv";
```

set output images path
``` C++
string downsampled_file_f  = "./Output/Downsampled_image_with_anti_aliasing.yuv";
string downsampled_file_nf = "./Output/Downsampled_image_without_anti_aliasing.yuv";
string output_file_f  = "./Output/Reconstructed_image_with_anti_aliasing.yuv";
string output_file_nf = "./Output/Reconstructed_image_without_anti_aliasing.yuv";
```

set input image's width and height
``` C++
int width  = 1920;
int height = 1080; 
```

set anti-aliasing level
``` C++
int sigma = 1;
```
and run program!

Then you can get below output:

+ 1/2 downsampled image with/without anti-aliasing

  ![image](https://user-images.githubusercontent.com/26856370/148549493-c4ac451e-7223-40ed-8614-d7eb75ca724d.png) 

+ reconstructed image with anti-aliasing

  ![image](https://user-images.githubusercontent.com/26856370/148549845-e3597e96-4ff5-4d71-beaf-49d38b0afec3.png)

+ reconstructed image without anti-aliasing

  ![image](https://user-images.githubusercontent.com/26856370/148549878-a5a7692b-08af-458e-adc5-eb170ff3ad46.png)

+ MSE of two reconstructed images

  ![image](https://user-images.githubusercontent.com/26856370/148550205-410ff7cb-8be3-4177-bc58-7549386d2a52.png)

