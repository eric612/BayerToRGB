#ifndef UTIL_H
#define UTIL_H

namespace RawScaler {

	typedef enum BayerType
	{
		RGGB = 0,
		GRBG = 1,
		BGGR = 2,
		GBRG = 3
	} ;
	float GetSum(float* arr, int num);
	void DivideArray(float* arr, int num, float divide);
	float MultiplyArrayAndSum(float* arr1, float* arr2, int num);
	float GetMean(float* arr, int num);
	float GetStd(float* arr, int num);
	float GetMedian(float* arr, int num);
	float GetMin(float* arr, int num);
	float GetMax(float* arr, int num);
	float distance(float x, float y, float i, float j);
	float pow_distance(float x, float y, float i, float j);
	double gaussian(float x, double sigma);
	unsigned char fixpoint_cut_uc(int value, int bit);
	char fixpoint_cut_c(int value, int bit);
	//void DownSample(unsigned char* in, unsigned char* out, int width, int height, BayerType type);
}
#endif