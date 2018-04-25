#include "util.h"
#include <math.h>
#include <string.h>
namespace RawScaler {
	unsigned char fixpoint_cut_uc(int value, int bit)
	{
		int max_val = 1 << bit;
		if (bit >= 8) {
			return 0;
		}
		if (value >= max_val) {
			value = max_val - 1;
		}
		unsigned char ans = value;
		return ans;

	}
	char fixpoint_cut_c(int value, int bit)
	{
		int max_val = 1 << (bit-1);
		if (bit-1 >= 8) {
			return 0;
		}
		if (value >= max_val) {
			value = max_val - 1;
		} 
		else if (value <= -max_val) {
			value = -(max_val - 1);
		}
		char ans = value;
		return ans;

	}
	float GetSum(float* arr, int num)
	{
		float sum = 0;
		for (int i = 0; i < num; i++) {
			sum += (arr[i]);
		}
		return sum;
	}
	void DivideArray(float* arr, int num, float divide)
	{
		float sum = 0;
		if (divide == 0) {
			return;
		}
		for (int i = 0; i < num; i++) {
			arr[i] = (arr[i] / divide);
		}
	}
	float MultiplyArrayAndSum(float* arr1, float* arr2, int num)
	{
		float* out = new float[num];
		float sum = 0;
		for (int i = 0; i < num; i++) {
			sum += (arr1[i] * arr2[i]);
		}
		delete[] out;
		return sum;
	}
	float GetMean(float* arr, int num)
	{
		float sum = 0;
		for (int i = 0; i < num; i++) {
			sum += (arr[i]);
		}
		sum = sum / (float)num;
		return sum;
	}
	float GetStd(float* arr, int num)
	{
		float mean = GetMean(arr, num);
		float sum = 0;
		for (int i = 0; i < num; i++) {
			sum += pow(arr[i] - mean, 2);
		}
		sum = sqrt(sum)/(float)num;
		return sum;
	}

	float GetMedian(float* arr, int num) {
		float temp;
		int i, j;
		// the following two loops sort the array x in ascending order
		for (i = 0; i<num - 1; i++) {
			for (j = i + 1; j<num; j++) {
				if (arr[j] < arr[i]) {
					// swap elements
					temp = arr[i];
					arr[i] = arr[j];
					arr[j] = temp;
				}
			}
		}

		if (num % 2 == 0) {
			// if there is an even number of elements, return mean of the two elements in the middle
			return((arr[num / 2] + arr[num / 2 - 1]) / 2.0);
		}
		else {
			// else return the element in the middle
			return arr[num / 2];
		}
	}
	float GetMin(float* arr, int num)
	{
		int idx = 0;
		float min_diff = 99999;
		for (int i = 0; i < num; i++) {
			float diff = fabs(arr[i]);
			if (diff < min_diff) {
				idx = i;
				min_diff = diff;
			}
		}
		return arr[idx];
	}
	float GetMax(float* arr, int num)
	{
		int idx = 0;
		float max_diff = 0;
		for (int i = 0; i < num; i++) {
			float diff = fabs(arr[i]);
			if (diff > max_diff) {
				idx = i;
				max_diff = diff;
			}
		}
		return arr[idx];
	}
	float distance(float x, float y, float i, float j) {
		return float(sqrt(pow(x - i, 2) + pow(y - j, 2)));
	}
	float pow_distance(float x, float y, float i, float j) {
		return float((fabs(x - i) * fabs(y - j)));
	}
	double gaussian(float x, double sigma) {
		return exp(-(pow(x, 2)) / (2 * pow(sigma, 2))) / (2 * 3.14159 * pow(sigma, 2));

	}

}