/*
 * Program's process:
 *  input image       -> downsampling -> upsampling -> reconstructed image
 *   ¤¤ anti-aliasing -> downsampling -> upsampling -> reconstructed image
 *
 *  You can check the effects of anti-aliasing by comparing output images and their MSE.
 *  This program is valid for YUV 4:2:0 8bit format images.
 */


#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

void ReadImage(int, int, string, char*, char*, char*);
void WriteImage(int, int, string, char*, char*, char*);
void Anti_aliasing(int, int, int, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*);
void CopyPadding(int, int, int, unsigned char**);
void DownSampling(int, int, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*);
void UpSampling(int, int, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*);
int GetMSE(int, int, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, unsigned char*);

int main()
{
	/* set input image's path */
	//string input_file	  = "./Input/PartyScene_832x480_yuv420_8bit_frame0.yuv";
	//string input_file	  = "./Input/RaceHorsesC_832x480_yuv420_8bit_frame0.yuv";
	//string input_file	  = "./Input/BasketballDrill_832x480_yuv420_8bit_frame0.yuv";
	string input_file = "./Input/Cactus_1920x1080_yuv420_8bit_frame200.yuv";
	/* set downsampled image's output path */
	string downsampled_file_f = "./Output/Downsampled_image_with_anti_aliasing.yuv";
	string downsampled_file_nf = "./Output/Downsampled_image_without_anti_aliasing.yuv";
	/* set reconstructed image's output path */
	string output_file_f = "./Output/Reconstructed_image_with_anti_aliasing.yuv";
	string output_file_nf = "./Output/Reconstructed_image_without_anti_aliasing.yuv";

	/* set input image's width & height */
	int width = 1920;
	int height = 1080;
	/* set gaussian filter(anti-aliasing filter)'s sigma */
	int sigma = 1;

	/* original data */
	unsigned char* Y = new unsigned char[width * height];
	unsigned char* U = new unsigned char[width * height / 4];
	unsigned char* V = new unsigned char[width * height / 4];
	/* anti-aliasing filtered data -> upsampled data with anti-aliasing */
	unsigned char* fY = new unsigned char[width * height];
	unsigned char* fU = new unsigned char[width * height / 4];
	unsigned char* fV = new unsigned char[width * height / 4];
	/* downsampled data */
	unsigned char* dY = new unsigned char[width * height / 4];
	unsigned char* dU = new unsigned char[width * height / 16];
	unsigned char* dV = new unsigned char[width * height / 16];
	/* upsampled data without anti-aliasing */
	unsigned char* rY = new unsigned char[width * height];
	unsigned char* rU = new unsigned char[width * height / 4];
	unsigned char* rV = new unsigned char[width * height / 4];

	/* Read original data */
	ReadImage(width, height, input_file, (char*)Y, (char*)U, (char*)V);

	/* Anti-aliasing filtering */
	Anti_aliasing(width, height, sigma, Y, U, V, fY, fU, fV);

	/* Downsampling & Upsampling original data */
	DownSampling(width, height, Y, U, V, dY, dU, dV);
	WriteImage(width / 2, height / 2, downsampled_file_nf, (char*)dY, (char*)dU, (char*)dV);
	UpSampling(width, height, dY, dU, dV, rY, rU, rV);
	WriteImage(width, height, output_file_nf, (char*)rY, (char*)rU, (char*)rV);

	/* Downsampling & Upsampling anti-aliasing filtered data */
	DownSampling(width, height, fY, fU, fV, dY, dU, dV);
	WriteImage(width / 2, height / 2, downsampled_file_f, (char*)dY, (char*)dU, (char*)dV);
	UpSampling(width, height, dY, dU, dV, fY, fU, fV);
	WriteImage(width, height, output_file_f, (char*)fY, (char*)fU, (char*)fV);

	/* Print MSE of reconstructed images */
	printf("MSE of reconstructed image without anti-aliasing : %d\n", GetMSE(width, height, Y, U, V, rY, rU, rV));
	printf("MSE of reconstructed image with anti-aliasing    : %d\n", GetMSE(width, height, Y, U, V, fY, fU, fV));

	delete[]Y, U, V;
	delete[]fY, fU, fV;
	delete[]dY, dU, dV;
	delete[]rY, rU, rV;
}


/* Read Y, Cu, Cv data of image */
void ReadImage(int width, int height,      /* image width & height */
			   string input_file,          /* file name to read    */
			   char* Y, char* U, char* V)  /* array to save data   */
{
	/* open image file */
	ifstream fin(input_file.c_str(), ios_base::in | ios_base::binary);
	if (!fin.is_open()) {                            
		printf("Error : Input file open fail"); 
		exit(1);                                      
	}       
	/* read image file */
	fin.read(Y, height * width);
	fin.read(U, height * width / 4);
	fin.read(V, height * width / 4);
	/* close image file */
	fin.close();
}
/* Write Y, Cu, Cv data of image */
void WriteImage(int width, int heigth,      /* image width & height */
				string output_file,         /* file name to write   */
				char* Y, char* U, char* V)  /* write data's array   */
{
	/* create image file */
	ofstream fout(output_file.c_str(), ios::out | ios::binary);
	if (!fout.is_open()) {                              
		printf("Error : Output file create fail");  
		exit(1);                                      
	}
	/* write image file */
	fout.write(Y, width * heigth);
	fout.write(U, width * heigth / 4);
	fout.write(V, width * heigth / 4);
	/* close image file */
	fout.close();
}

/* Get anti-aliasing filtered image data from input image data using 1D gaussian filter */
void Anti_aliasing(int width, int height,   /* image width & height    */
	               int sigma,               /* sigma value of gaussian */
				   unsigned char* Y,        /* input image Y data      */
				   unsigned char* Cu,       /* input image U data      */
				   unsigned char* Cv,       /* input image V data      */
				   unsigned char* fY,       /* filtered image Y data   */ 
				   unsigned char* fCu,      /* filtered image U data   */ 
				   unsigned char* fCv)      /* filtered image V data   */ 
{
	/* set kernel & padding size */
	int kernel_size = 6 * sigma + 1;
	int padding_size = 3 * sigma;
	/* array to save 1D gaussian filter */
	unsigned int* kernel = new unsigned int[kernel_size];

	int sum = 0, idx = 0;
	double gaussian;
	/* get edge value of 1D gaussian filter */
	double edge = 1 / (sqrt(2 * M_PI) * sigma) * pow(M_E, (-1 * pow(padding_size, 2) / (2 * pow(sigma, 2))));
	/* get 1D gaussian filter */ 
	// loop half of kernel size (gaussian is symmetric)
	for (int i = padding_size; i >= 0; i--)
	{
		// get 1D gaussian function value
		gaussian = 1 / (sqrt(2 * M_PI) * sigma) * pow(M_E, (-1 * pow(i, 2) / (2 * pow(sigma, 2))));
		// normalize to integer
		kernel[idx] = floor(gaussian / edge);
		// save to symmetric position
		kernel[kernel_size - 1 - idx] = kernel[idx];
		// get normalization weight
		sum += (2 * kernel[idx]);
		idx++;
	}
	sum -= kernel[kernel_size / 2];

	/* 2D array to save image data with padding */
	unsigned char** Y_padding  = new unsigned char* [height + 2 * padding_size];
	unsigned char** Cu_padding = new unsigned char* [height / 2 + 2 * padding_size];
	unsigned char** Cv_padding = new unsigned char* [height / 2 + 2 * padding_size];
	for (int i = 0; i < height + 2 * padding_size; i++)
	{
		Y_padding[i] = new unsigned char[width + 2 * padding_size];
		if (i < height / 2 + 2 * padding_size) {
			Cu_padding[i] = new unsigned char[width / 2 + 2 * padding_size];
			Cv_padding[i] = new unsigned char[width / 2 + 2 * padding_size];
		}
	}
	/* copy image data from 1D array to 2D array */
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			Y_padding[i + padding_size][j + padding_size] = Y[i * width + j];
			if (i < height / 2 && j < width / 2) {
				Cu_padding[i + padding_size][j + padding_size] = Cu[i * width / 2 + j];
				Cv_padding[i + padding_size][j + padding_size] = Cv[i * width / 2 + j];
			}
		}
	}
	/* fill padding data */
	CopyPadding(width,     height,     padding_size, Y_padding );
	CopyPadding(width / 2, height / 2, padding_size, Cu_padding);
	CopyPadding(width / 2, height / 2, padding_size, Cv_padding);
	
	/* 2D array to save intermediate image data after first 1D filtering */
	unsigned char** Y_temp  = new unsigned char* [height];
	unsigned char** Cu_temp = new unsigned char* [height / 2];
	unsigned char** Cv_temp = new unsigned char* [height / 2];
	for (int i = 0; i < height; i++)
	{
		Y_temp[i] = new unsigned char[width + 2 * padding_size];
		if (i < height / 2) {
			Cu_temp[i] = new unsigned char[width / 2 + 2 * padding_size];
			Cv_temp[i] = new unsigned char[width / 2 + 2 * padding_size];
		}
	}

	/* first 1D filtering in column direction */
	for (int i = padding_size; i < height + padding_size; i++)
	{
		for (int j = 0; j < width + 2 * padding_size; j++)
		{
			// convolution of 1D gaussian filter and Y data
			int value = 0;
			for (int k = 0; k < kernel_size; k++)
				value += (kernel[k] * Y_padding[i - padding_size + k][j]);
			Y_temp[i - padding_size][j] = value / sum;

			// convolution of 1D gaussian filter and U data 
			// convolution of 1D gaussian filter and V data 
			if (i < height / 2 + padding_size && j < width / 2 + 2 * padding_size) {
				int value1 = 0, value2 = 0;
				for (int k = 0; k < kernel_size; k++)
				{
					value1 += (kernel[k] * Cu_padding[i - padding_size + k][j]);
					value2 += (kernel[k] * Cv_padding[i - padding_size + k][j]);
				}
				Cu_temp[i - padding_size][j] = value1 / sum;
				Cv_temp[i - padding_size][j] = value2 / sum;
			}
		}
	}
	/* second 1D filtering in row direction */
	for (int i = 0; i < height; i++)
	{
		for (int j = padding_size; j < width + padding_size; j++)
		{
			// convolution of 1D gaussian filter and Y data
			int value = 0;
			for (int k = 0; k < kernel_size; k++)
				value += (kernel[k] * Y_temp[i][j - padding_size + k]);
			fY[i * width + j - padding_size] = value / sum;;

			// convolution of 1D gaussian filter and U data 
			// convolution of 1D gaussian filter and V data 
			if (i < height / 2 && j < width / 2 + padding_size) {
				int value1 = 0, value2 = 0;
				for (int k = 0; k < kernel_size; k++)
				{
					value1 += (kernel[k] * Cu_temp[i][j - padding_size + k]);
					value2 += (kernel[k] * Cv_temp[i][j - padding_size + k]);
				}
				fCu[i * width / 2 + j - padding_size] = value1 / sum;
				fCv[i * width / 2 + j - padding_size] = value2 / sum;
			}
		}
	}

	for (int i = 0; i < height + 2 * padding_size; i++)
		delete[]Y_padding[i];
	for (int i = 0; i < height / 2 + 2 * padding_size; i++)
	{
		delete[]Cu_padding[i];
		delete[]Cv_padding[i];
	}
	delete[]kernel, Y_padding, Cu_padding, Cv_padding;
}
/* Get image data with replicate(copy) padding */
void CopyPadding(int width, int height,  /* image width & height  */
				 int padding_size,       /* padding size          */
				 unsigned char** input)  /* image data's 2D array */
{
	/* vertical padding */
	for (int j = padding_size; j < width + padding_size; j++)
	{
		for (int i = 0; i < padding_size; i++)
		{
			// upward
			input[i][j] = input[padding_size][j];
			// downward
			input[i + height + padding_size][j] = input[height + padding_size - 1][j];
		}
	}
	/* horiontal(+ diagonal) padding */
	for (int i = 0; i < height + 2 * padding_size; i++)
	{
		for (int j = 0; j < padding_size; j++)
		{
			// left direction 
			input[i][j] = input[i][padding_size];
			// right direction
			input[i][j + width + padding_size] = input[i][width + padding_size - 1];
		}
	}
}

/* Get downsampled image data from original image data */
void DownSampling(int width, int height,  /* image width & height     */
				  unsigned char* Y,       /* input image Y data       */
				  unsigned char* Cu,      /* input image U data       */
				  unsigned char* Cv,      /* input image V data       */
				  unsigned char* dY,      /* downsampled image Y data */
				  unsigned char* dCu,     /* downsampled image U data */
				  unsigned char* dCv)     /* downsampled image V data */
{
	/* downsampling to 1/2 width and 1/2 height */
	for (int i = 0; i < height / 2; i++)
	{
		for (int j = 0; j < width / 2; j++)
		{
			// save last value of each 2*2 window in input Y data
			dY[i * width / 2 + j] = Y[(2 * i + 1) * width + (2 * j + 1)];

			// save last value of each 2*2 window in input U data
			// save last value of each 2*2 window in input V data
			if (i < height / 4 && j < width / 4) {
				dCu[i * width / 4 + j] = Cu[(2 * i + 1) * width / 2 + (2 * j + 1)];
				dCv[i * width / 4 + j] = Cv[(2 * i + 1) * width / 2 + (2 * j + 1)];
			}
		}
	}
}

/* Get upsamplied image data from downsampled image data */
void UpSampling(int width, int height,  /* image width & height     */
				unsigned char* dY,      /* downsampled image Y data */
				unsigned char* dCu,     /* downsampled image U data */
				unsigned char* dCv,     /* downsampled image V data */
				unsigned char* Y,       /* upsampled image Y data   */
				unsigned char* Cu,      /* upsampled image U data   */
				unsigned char* Cv)      /* upsampled image V data   */ 
{
	/* Get integer-pel sample from downsampled image data */
	for (int i = 0; i < height / 2; i++)
	{
		for (int j = 0; j < width / 2; j++)
		{
			// save to last value of each 2*2 window in input Y data
			Y[(2 * i + 1) * width + (2 * j + 1)] = dY[i * width / 2 + j];

			// save to last value of each 2*2 window in input U data
			// save to last value of each 2*2 window in input V data
			if (i < height / 4 && j < width / 4) {
				Cu[(2 * i + 1) * width / 2 + (2 * j + 1)] = dCu[i * width / 4 + j]  ;
				Cv[(2 * i + 1) * width / 2 + (2 * j + 1)] = dCv[i * width / 4 + j];
			}
		}
	}

	/* Get horizaontal direction's half-pel sample */
	for (int i = 1; i < height; i += 2)
	{
		for (int j = 0; j < width; j += 2)
		{
			// Get  6 Y data used in interpolation (consider padding)
			int idx = 0;
			unsigned char temp[6];
			for (int k = -5; k <= 5; k += 2)
			{
				if (j + k < 0)
					temp[idx++] = Y[i * width + 1];
				else if(j + k >= width)
					temp[idx++] = Y[i * width + width - 1];
				else 
					temp[idx++] = Y[i * width + j + k];
			}
			// Do Y data horizaontal interpolation using 6-tap DCT-IF
			Y[i * width + j] = ((temp[0] + temp[5]) * 11 - (temp[1] + temp[4]) * 43 + (temp[2] + temp[3]) * 160) / 256;

			// Get 6 U data used in interpolation (consider padding)
			// Get 6 V data used in interpolation (consider padding)
			if (i < height / 2 && j < width / 2) {
				int idx1 = 0, idx2 = 0;
				unsigned char temp1[6], temp2[6];
				for (int k = -5; k <= 5; k += 2)
				{
					if (j + k < 0) {
						temp1[idx1++] = Cu[i * width / 2 + 1];
						temp2[idx2++] = Cv[i * width / 2 + 1];
					}
					else if (j + k >= width) {
						temp1[idx1++] = Cu[i * width / 2 + width / 2 - 1];
						temp2[idx2++] = Cv[i * width / 2 + width / 2 - 1];
					}
					else {
						temp1[idx1++] = Cu[i * width / 2 + j + k];
						temp2[idx2++] = Cv[i * width / 2 + j + k];
					}
				}
				// Do U data horizaontal interpolation using 6-tap DCT-IF
				Cu[i * width / 2 + j] = ((temp1[0] + temp1[5]) * 11 - (temp1[1] + temp1[4]) * 43 + (temp1[2] + temp1[3]) * 160) / 256;
				// Do V data horizaontal interpolation using 6-tap DCT-IF
				Cv[i * width / 2 + j] = ((temp2[0] + temp2[5]) * 11 - (temp2[1] + temp2[4]) * 43 + (temp2[2] + temp2[3]) * 160) / 256;
			}
		}
	}

	/* Get vertical direction's half-pel sample */
	for (int i = 0; i < height; i += 2)
	{
		for (int j = 1; j < width; j += 2)
		{
			// Get 6 Y data used in interpolation (consider padding)
			int idx = 0;
			unsigned char temp[6];
			for (int k = -5; k <= 5; k += 2)
			{
				if (i + k < 0)
					temp[idx++] = Y[width + j];
				else if (i + k >= height)
					temp[idx++] = Y[(height - 1) * width + j];
				else
					temp[idx++] = Y[(i + k) * width + j];
			}
			// Do Y data vertical interpolation using 6-tap DCT-IF
			Y[i * width + j] = ((temp[0] + temp[5]) * 11 - (temp[1] + temp[4]) * 43 + (temp[2] + temp[3]) * 160) / 256;

			// Get 6 U data used in interpolation (consider padding)
			// Get 6 V data used in interpolation (consider padding)
			if (i < height / 2 && j < width / 2) {
				int idx1 = 0, idx2 = 0;
				unsigned char temp1[6], temp2[6];
				for (int k = -5; k <= 5; k += 2)
				{
					if (i + k < 0) {
						temp1[idx1++] = Cu[width / 2 + j];
						temp2[idx2++] = Cv[width / 2 + j];
					}
					else if (i + k >= width) {
						temp1[idx1++] = Cu[(height / 2 - 1) * width / 2 + j];
						temp2[idx2++] = Cv[(height / 2 - 1) * width / 2 + j];
					}
					else {
						temp1[idx1++] = Cu[(i + k) * width / 2 + j];
						temp2[idx2++] = Cv[(i + k) * width / 2 + j];
					}
				}
				// Do U data vertical interpolation using 6-tap DCT-IF
				Cu[i * width / 2 + j] = ((temp1[0] + temp1[5]) * 11 - (temp1[1] + temp1[4]) * 43 + (temp1[2] + temp1[3]) * 160) / 256;
				// Do V data vertical interpolation using 6-tap DCT-IF
				Cv[i * width / 2 + j] = ((temp2[0] + temp2[5]) * 11 - (temp2[1] + temp2[4]) * 43 + (temp2[2] + temp2[3]) * 160) / 256;
			}
		}
	}

	/* Get diagonal direction's half-pel sample 
	  (horizontal interpolation using vertical half-pel sampel) */
	for (int i = 0; i < height; i += 2)
	{
		for (int j = 0; j < width; j += 2)
		{
			// Get 6 Y data used in interpolation (consider padding)
			int idx = 0;
			unsigned char temp[6];
			for (int k = -5; k <= 5; k += 2)
			{
				if (j + k < 0)
					temp[idx++] = Y[i * width + 1];
				else if (j + k >= width)
					temp[idx++] = Y[i * width + width - 1];
				else
					temp[idx++] = Y[i * width + j + k];
			}
			// Do Y data horizontal interpolation using 6-tap DCT-IF
			Y[i * width + j] = ((temp[0] + temp[5]) * 11 - (temp[1] + temp[4]) * 43 + (temp[2] + temp[3]) * 160) / 256;

			// Get 6 U data used in interpolation (consider padding)
			// Get 6 V data used in interpolation (consider padding)
			if (i < height / 2 && j < width / 2) {
				int idx1 = 0, idx2 = 0;
				unsigned char temp1[6], temp2[6];
				for (int k = -5; k <= 5; k += 2)
				{
					if (j + k < 0) {
						temp1[idx1++] = Cu[i * width / 2 + 1];
						temp2[idx2++] = Cv[i * width / 2 + 1];
					}
					else if (j + k >= width) {
						temp1[idx1++] = Cu[i * width / 2 + width / 2 - 1];
						temp2[idx2++] = Cv[i * width / 2 + width / 2 - 1];
					}
					else {
						temp1[idx1++] = Cu[i * width / 2 + j + k];
						temp2[idx2++] = Cv[i * width / 2 + j + k];
					}
				}
				// Do U data horizontal interpolation using 6-tap DCT-IF
				Cu[i * width / 2 + j] = ((temp1[0] + temp1[5]) * 11 - (temp1[1] + temp1[4]) * 43 + (temp1[2] + temp1[3]) * 160) / 256;
				// Do V data horizontal interpolation using 6-tap DCT-IF
				Cv[i * width / 2 + j] = ((temp2[0] + temp2[5]) * 11 - (temp2[1] + temp2[4]) * 43 + (temp2[2] + temp2[3]) * 160) / 256;
			}
		}
	}
}

/* Get Mean Square Error from two image data */
int GetMSE(int width, int heigth,  /* image width & height */
		   unsigned char* Y,       /* image1 Y data        */
		   unsigned char* U,       /* image1 U data        */
		   unsigned char* V,       /* image1 V data        */
		   unsigned char* rY,      /* image2 Y data        */
		   unsigned char* rU,      /* image2 U data        */
		   unsigned char* rV)      /* image2 V data        */
{
	int MSE_Y = 0, MSE_UV = 0;
	for (int i = 0; i < heigth; i++)
	{
		for (int j = 0; j < width; j++)
		{
			// Calculate Y data's SE
			MSE_Y += (pow(Y[i * width + j] - rY[i * width + j], 2));

			// Calculate U data's SE
			// Calculate V data's SE
			if (i < heigth / 2 && j < width / 2)
			{
				MSE_UV += (pow(U[i * width / 2 + j] - rU[i * width / 2 + j], 2));
				MSE_UV += (pow(V[i * width / 2 + j] - rV[i * width / 2 + j], 2));
			}
		}
	}
	return (MSE_Y / (width * heigth)) + (MSE_UV / (width * heigth / 4));
}
