// A simple program that computes the square root of a number
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef NoCmake
#include "RawScaler.config.h"
#endif

#ifndef OpenCV
//#include "BasicBitmap.h"
#include "bitmap_image.hpp"
#else
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#ifdef openmp
#include <omp.h>
#endif
#endif // !OpenCV
#include "rawscaler.h"
#include <math.h>
#include <time.h>   /* 時間相關函數 */
using namespace RawScaler;  

BayerType mBayerType = BGGR; //GRBG,BGGR,RGGB,GBRG
#define BOUND(a,min_val,max_val)           ( (a < min_val) ? min_val : (a >= max_val) ? (max_val) : a )

#ifdef OpenCV
cv::Mat RGB2RAW(cv::Mat src, BayerType type = BGGR)
{
	
	int w = src.cols;
	int h = src.rows;
	cv::Mat raw_img(h,w, CV_8UC1);
#pragma omp parallel for
	for (int i = 0; i < h; i++) {
		const uchar* ptr1 = src.ptr<uchar>(i);
		uchar* ptr2 = raw_img.ptr<uchar>(i);
		int img_index1 = 0;
		int img_index2 = 0;
		int offset = i % 2;
		for (int j = 0; j < w; j++)
		{
			int ch = (j) % 2;
			if (!offset) { 
				switch (type)
				{
				case RGGB:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 2];
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 1]; // G channel
					}
					break;
				case GRBG:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 1];
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 2]; // G channel
					}
					break;
				case BGGR:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 0];
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 1]; // G channel
					}
					break;
				case GBRG:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 1];
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 0]; // G channel
					}
					break;
				default:
					printf("error bayer type\n");
					//return raw_img;
				}

			}
			else {
				switch (type)
				{
				case RGGB:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 1]; // G channel
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 0];
					}
					break;
				case GRBG:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 0]; // G channel
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 1];
					}
					break;
				case BGGR:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 1]; // G channel
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 2];
					}
					break;
				case GBRG:
					if (!ch) {
						ptr2[img_index2] = ptr1[img_index1 + 2]; // G channel
					}
					else {
						ptr2[img_index2] = ptr1[img_index1 + 1];
					}
					break;
				default:
					printf("error bayer type\n");
					//return raw_img;
				}

			}
			//float l = pow(ptr2[img_index2]/255.0,  2.2f);
			//ptr2[img_index2] = l * 255;
			img_index2++;
			img_index1 += 3;
		}
	}
	return raw_img;
}
template <typename Dtype>
void Mat2Char(cv::Mat src, Dtype* out=NULL, int channel = 1,int multiply = 1,bool randval = false)
{
	uchar depth = src.type() & CV_MAT_DEPTH_MASK;
	if (out == NULL) {
		printf("Mat2Char,output image memory was NULL\n");
		return;
	}
	else if (depth != CV_8UC1) {
		printf("Mat2Char,source mat was not CV_8UC1\n");
		return;
	}
	int i, j;
	int height = src.rows;
	int width = src.cols;
	srand(time(NULL));
	int max_val = 256 * multiply - 1;
	float rand_val = ((rand() % 25) / 100.);
	int rand_sign = rand() % 2;
	if (rand_sign == 1) {
		rand_val = -rand_val;
	}
//#pragma omp parallel for
	for (i = 0; i < height; i++) {
		const uchar* ptr = src.ptr<uchar>(i);
		int img_index = 0;
		for (j = 0; j < width*channel; j++) {
			if (randval) {
				int rand_val2 = ((rand() % 5))*multiply;
				int rand_sign2 = rand() % 2;
				if (rand_sign2 == 1) {
					rand_val2 = -rand_val2;
				}
				out[i*width*channel + img_index] = BOUND(ptr[img_index] * (multiply*(1+rand_val) )+ rand_val2, 0, max_val);
			}	                        
			else {
				out[i*width*channel + img_index] = ptr[img_index]* multiply;
			}
			img_index++;
		}
	}
}
template <typename Dtype>
void Char2Mat(Dtype* src, cv::Mat out,int channel = 1, int multiply = 1)
{
	uchar depth = out.type() & CV_MAT_DEPTH_MASK;
	if (src == NULL) {
		printf("Mat2Char,output image memory was NULL\n");
		return;
	}
	else if (depth != CV_8UC1) {
		printf("Mat2Char,source mat was not CV_8UC1\n");
		return;
	}
	int i, j;
	int height = out.rows;
	int width = out.cols;
//#pragma omp parallel for
	for (i = 0; i < height; i++) {
		uchar* ptr = out.ptr<uchar>(i);
		int img_index = 0;
		for (j = 0; j < width*channel; j++) {
			ptr[img_index] = src[i*width*channel + img_index]/ multiply;
			img_index++;
		}
	}
}
cv::Mat ConvertBayer2BGR(cv::Mat bayer)
{
	cv::Mat BGR;// (bayer.rows, bayer.cols, CV_8UC3);

	switch (mBayerType)
	{
	case RGGB:
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerBG2BGR);//COLOR_BayerBG2BGR, COLOR_BayerBG2GRAY
		break;
	case GRBG:
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerGB2BGR);
		break;
	case BGGR:	
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerRG2BGR);
		break;
	case GBRG:
		
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerGR2BGR);
		break;
	default:
		printf("error bayer type\n");
	}
	
	//Char2Mat(scaler_img2, BGR);
	return BGR;
}
cv::Mat ConvertBayer2GRAY(cv::Mat bayer)
{
	cv::Mat BGR;
	
	switch (mBayerType)
	{
	case RGGB:
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerBG2GRAY);
		break;
	case GRBG:
		cv::cvtColor(bayer, BGR, cv::COLOR_BayerGB2GRAY);
		break;
	case BGGR:

		cv::cvtColor(bayer, BGR, cv::COLOR_BayerRG2GRAY);	
		break;
	case GBRG:

		cv::cvtColor(bayer, BGR, cv::COLOR_BayerGR2GRAY);
		break;
	default:
		printf("error bayer type\n");
	}

	return BGR;
}
cv::Mat MakeRawBorder(cv::Mat img1,int pad_width)
{
	int i, j;
	double eqm = 0;
	int height = img1.rows;
	int width = img1.cols;
	cv::Mat img2(height + pad_width *2, width+pad_width*2, CV_8UC1);
//#pragma omp parallel for
	for (i = 0; i < height; i++) {
		const uchar* ptr1 = img1.ptr<uchar>(i);
		uchar* ptr2 = img2.ptr<uchar>(i+pad_width);
		int img_index1 = 0;
		int img_index2 = pad_width;
		for (j = 0; j < width; j++) {
			
			if (j == 0) {
				for (int k = 0; k < pad_width + 1; k++) {
					if (k % 2 == 0) {
						ptr2[img_index2 - k] = ptr1[img_index1];
					}
					else {
						ptr2[img_index2 - k] = ptr1[img_index1 + 1];
					}
				}
			}
			else if(j==width-1) {
				for (int k = 0; k < pad_width + 1; k++) {
					
					if (k % 2 == 0) {
						ptr2[img_index2 + k] = ptr1[img_index1];
					}
					else {
						ptr2[img_index2 + k] = ptr1[img_index1 - 1];
					}
				}
			}
			else {
				ptr2[img_index2] = ptr1[img_index1];
			}
			img_index1++;
			img_index2++;
		}
	}
	for (i = 0; i < pad_width; i++) {
		const uchar* ptr1;
		if (i % 2 == 0) {
			ptr1 = img2.ptr<uchar>(pad_width+1);
		}
		else {
			ptr1 = img2.ptr<uchar>(pad_width);
		}
		
		uchar* ptr2 = img2.ptr<uchar>(pad_width - 1 -i);
		int img_index1 = 0;
		int img_index2 = 0;
		for (j = 0; j < width + pad_width*2; j++) {
			ptr2[img_index2] = ptr1[img_index1];
			img_index1++;
			img_index2++;
		}
	}
	for (i = 0; i < pad_width; i++) {
		const uchar* ptr1;
		if (i % 2 == 0) {
			ptr1 = img2.ptr<uchar>(height + pad_width - 2);
		}
		else {
			ptr1 = img2.ptr<uchar>(height + pad_width - 1);
		}

		uchar* ptr2 = img2.ptr<uchar>(pad_width + height + i);
		int img_index1 = 0;
		int img_index2 = 0;
		for (j = 0; j < width + pad_width * 2; j++) {
			ptr2[img_index2] = ptr1[img_index1];
			img_index1++;
			img_index2++;
		}
	}
	//eqm /= height * width;
	return img2;
}
void Bayer2RGB(unsigned char* img,int in_width,int in_height,int out_width,int out_height)
{
	char dir[128];
	sprintf(dir, "%s", OutputDir);
	RawProcessor scaler;
	cv::Mat dst;
	cv::Mat src(in_height, in_width, CV_8UC1);
	Char2Mat(img, src, 1);
	//cv::imwrite("../../out/raw_scaler.png", src);
	int pad = 6;
	int pad_w = pad * 2;
	dst = MakeRawBorder(src, pad);
	char filename[128];
	sprintf(filename, "%spad.bmp", dir);

	cv::imwrite(filename, dst);
	unsigned char* dw_img4 = new unsigned char[(in_width + pad_w)*(in_height + pad_w)];
	unsigned char* dw_img5 = new unsigned char[(out_width + pad_w)*(out_height + pad_w) * 3];
	Mat2Char(dst, dw_img4, 1);

	scaler.BayerToRGB(dw_img4, dw_img5, in_width + pad_w, in_height + pad_w, out_width + pad_w, out_height + pad_w, mBayerType);


	cv::Mat test(out_height + pad_w, out_width + pad_w, CV_8UC3);

	Char2Mat(dw_img5, test, 3);

	cv::Rect myROI(pad, pad, out_width, out_height);
	// Crop the full image to that image contained by the rectangle myROI
	// Note that this doesn't copy the data
	cv::Mat croppedImage = test(myROI);

	sprintf(filename, "%sresult.png", dir);

	cv::imwrite(filename, croppedImage);

	cv::imshow("demo", croppedImage);
	src.release();
	test.release();
	delete[] dw_img4;
	delete[] dw_img5;
}
int main(int argc, char *argv[])
{
	int res = -1;
	if (argc > 1) {
		res = strcmp(argv[1], "test");
	}
	char dir[128];
	sprintf(dir, "%s", OutputDir);

#ifdef openmp
	printf("\nuse openmp\n");
#endif
	cv::Mat img,res_img,truth;

	int user_w = 0;
	int user_h = 0;
	if (argc >= 2)
	{
		img = cv::imread(argv[1]);
		if (argc >= 4) {
			user_w = atoi(argv[2]);
			user_h = atoi(argv[3]);
		}
	}

	if (img.data == NULL) {
		img = cv::imread(ImagePath);
	}
	RawProcessor scaler;
	cv::Mat raw_img = RGB2RAW(img, mBayerType);

	float ratio = 1.;
	int in_height = raw_img.rows ;
	int in_width = raw_img.cols ;
	
	int out_height;// = raw_img.rows / ratio;
	int out_width;// = raw_img.cols / ratio;
	if (user_w <= in_width && user_w > 32 && user_h <= in_height && user_h > 32) {
		out_width = user_w;
		out_height = user_h;
	}
	else {
		out_height = raw_img.rows / ratio;
		out_width = raw_img.cols / ratio;

		out_width = out_width + (out_width % 2);
		out_height = out_height + (out_height % 2);
	}
	unsigned char* scaler_img = new unsigned char[raw_img.cols*raw_img.rows];
	cv::Mat raw_img2 = RGB2RAW(img, mBayerType);
	Mat2Char(raw_img2, scaler_img);

	Bayer2RGB(scaler_img, in_width, in_height, out_width, out_height);
	cv::waitKey(3000);
	
	delete[] scaler_img;
	img.release();

	raw_img.release();
	//res_img.release();
	truth.release();
	return 0;
}
#else
unsigned char* RGB2RAW(bitmap_image src, BayerType type = BGGR)
{
	int w = src.get_width();
	int h = src.get_height();
	//#pragma omp parallel for
	unsigned char* int_img = new unsigned char[w*h];
	for (int i = 0; i < h; i++) {
		int img_index1 = 0;
		int img_index2 = i*w;
		int offset = i % 2;
		for (int j = 0; j < w; j++)
		{
			int ch = (j) % 2;
			//IUINT32 clr = src->GetPixel(j, i);
			rgb_t t = src.get_pixel(j, i);
			//int red = clr & 0xff;
			//int green = (clr >> 8) & 0xff;
			//int blue = (clr >> 16) & 0xff;
			if (!offset) {
				switch (type)
				{
				case RGGB:
					if (!ch) {
						int_img[img_index2] = t.red;
					}
					else {
						int_img[img_index2] = t.green;
					}
					break;
				case GRBG:
					if (!ch) {
						int_img[img_index2] = t.green;
					}
					else {
						int_img[img_index2] = t.red;
					}
					break;
				case BGGR:
					if (!ch) {
						int_img[img_index2] = t.blue;
					}
					else {
						int_img[img_index2] = t.green;
					}
					break;
				case GBRG:
					if (!ch) {
						int_img[img_index2] = t.green;
					}
					else {
						int_img[img_index2] = t.blue;
					}
					break;
				default:
					printf("error bayer type\n");
					//return raw_img;
				}

			}
			else {
				switch (mBayerType)
				{
				case RGGB:
					if (!ch) {
						int_img[img_index2] = t.green;
					}
					else {
						int_img[img_index2] = t.blue;
					}
					break;
				case GRBG:
					if (!ch) {
						int_img[img_index2] = t.blue;
					}
					else {
						int_img[img_index2] = t.green;
					}
					break;
				case BGGR:
					if (!ch) {
						int_img[img_index2] = t.green;
					}
					else {
						int_img[img_index2] = t.red;
					}
					break;
				case GBRG:
					if (!ch) {
						int_img[img_index2] = t.red;
					}
					else {
						int_img[img_index2] = t.green;
					}
					break;
				default:
					printf("error bayer type\n");
					//return raw_img;
				}

			}
			//float l = pow(ptr2[img_index2]/255.0,  2.2f);
			//ptr2[img_index2] = l * 255;
			img_index2++;
			img_index1 += 3;
		}
	}
	return int_img;
}

int main(int argc, char *argv[])
{
	int res = -1;
	if (argc > 1) {
		res = strcmp(argv[1], "test");
	}
	char dir[128];
	sprintf(dir, "%s", OutputDir);

	res = -1;
	if (argc > 1) {
		res = strcmp(argv[1], "demo");
	}
	int index = 1;
	int user_w = 0;
	int user_h = 0;

	if (argc >= 2 && res == 0) {			
		if (argc >= 4) {
			user_w = atoi(argv[2]);
			user_h = atoi(argv[3]);
		}
	}
	bitmap_image image= bitmap_image(ImagePath);
	int in_width = image.get_width();
	int in_height = image.get_height();
	unsigned char* int_img = RGB2RAW(image, mBayerType);
	unsigned char* dw_img2 = new unsigned char[in_width*in_height];
	for (int i = 0; i < in_width*in_height; i++) {
		dw_img2[i] = int_img[i];
	}
	float ratio = 3.5;
	int out_height;// = raw_img.rows / ratio;
	int out_width;// = raw_img.cols / ratio;
	if (user_w <= in_width && user_w > 32 && user_h <= in_height && user_h > 32) {
		out_width = user_w;
		out_height = user_h;
	}
	else {
		out_height = in_height / ratio;
		out_width = in_width / ratio;

		out_width = out_width + (out_width % 2);
		out_height = out_height + (out_height % 2);
	}
	RawProcessor scaler;
	unsigned char* dw_img4 = new unsigned char[out_width*out_height*3];
	scaler.BayerToRGB(dw_img2, dw_img4, in_width, in_height, out_width, out_height, mBayerType);
	bitmap_image out_image(out_width, out_height);
	int w = out_width;
	for (int i = 0; i < out_height; i++) {
		int offset = i % 2;
		for (int j = 0; j < out_width; j++) {
			int val1 = dw_img4[i*out_width*3 + j*3];
			int val2 = dw_img4[i*out_width * 3 + j * 3 + 1];
			int val3 = dw_img4[i*out_width * 3 + j * 3 + 2];
			out_image.set_pixel(j, i, val3, val2, val1);
		}
	}
	sprintf(dir, "%s", OutputDir);
	out_image.save_image(strcat(dir,"result_raw.bmp"));
	delete[] int_img;
	delete[] dw_img2;
	//delete[] dw_img3;
	delete[] dw_img4;
	out_image.clear();
	image.clear();

	return 0;
}
#endif