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
	void ReSize(int iwidth,int iheight,int newwidth, int newheight,unsigned char* data, unsigned char* out)
	{

		int x, y;
		int y1, x1;
		double w1, h1;
		double x1_pos;
		double y1_pos;
		double dx, dy;
		double x1_offset;
		double y1_offset;
		int w_line_pixel;
		int i_pos;
		int o_pos;
		int m_iwidth = iwidth;
		int m_iheight = iheight;
		int i_width = m_iwidth;
		int i_height = m_iheight;
		int o_width = newwidth;
		int o_height = newheight;
		int m_ipaddedwidth = newwidth * 3;// +(newwidth & 3);
		unsigned char* i_data = new unsigned char[(i_width * 3 )*i_height];
		memcpy(i_data, data, (i_width * 3 )*i_height);
		//Delete();

		int padding = 0;
		unsigned char* m_pbitdata = new unsigned char[o_width*o_height * 3 ];

		m_iwidth = o_width;
		m_iheight = o_height;



		unsigned char* o_data = out;
		//initial value
		w_line_pixel = i_width * 3 ;
		i_pos = 0;
		o_pos = 0;


		//input image one pixel's w=1 and h=1
		//output image one pixel's w1 and h1
		//w*i_width = w1*o_width;
		w1 = (double)i_width / (double)o_width;
		h1 = (double)i_height / (double)o_height;

		x1_offset = w1 / 2.0;
		y1_offset = h1 / 2.0;

		double Weight[16];
		int k;
		// 00, 01, 02, 03
		// 04, 05, 06, 07
		// 08, 09, 10, 11
		// 12, 13, 14, 15

		Weight[0] = -5;
		Weight[1] = -3;
		Weight[2] = -3;
		Weight[3] = -5;

		Weight[4] = -3;
		Weight[7] = -3;

		double v = 60. / 67.;
		Weight[5] = v;
		Weight[6] = v;
		Weight[9] = v;
		Weight[10] = v;

		Weight[8] = -3;
		Weight[11] = -3;

		Weight[12] = -5;
		Weight[13] = -3;
		Weight[14] = -3;
		Weight[15] = -5;

		for (y1 = 0; y1<o_height; y1++)
		{

			y1_pos = y1*h1 + y1_offset;
			y = int(y1_pos - 0.5);
			dy = y1_pos - y - 0.5;

			if (dy <= 0)
				dy = 0.0000001;
			if (y >= i_height)
				y = i_height - 1;

			for (x1 = 0; x1<o_width; x1++)
			{

				x1_pos = x1*w1 + x1_offset;
				x = int(x1_pos - 0.5);
				dx = x1_pos - x - 0.5;
				if (dx <= 0)
					dx = 0.0000001;

				if (x >= i_width)
					x = i_width - 1;

				double R, G, B;
				double ctt;
				double weight_val;
				double area;
				ctt = 0;
				R = G = B = 0;

				k = 0;
				//start filter operation
				for (int ty = -1; ty <= 2; ty++) {
					if ((ty + y)< 0 || (ty + y) >= i_height) {
						k += 4;
						continue;
					}
					i_pos = (y + ty)*w_line_pixel + (x - 1) * 3;

					for (int tx = -1; tx <= 2; tx++) {
						if ((tx + x)< 0 || (tx + x) >= i_width) {
							i_pos += 3;
							k++;
							continue;
						}
						//k=(ty+1)*4+(tx+1);
						area = (dx - tx)*(dy - ty);
						area = (area>0 ? area : -area);
						weight_val = 1.0 / (Weight[k++] * area);
						B += (weight_val*(double)i_data[i_pos++]);
						G += (weight_val*(double)i_data[i_pos++]);
						R += (weight_val*(double)i_data[i_pos++]);
						ctt += (weight_val);
					}
				}
				//end filter operation
				R /= ctt;
				R = (R>255 ? 255 : R);
				R = (R<0 ? 0 : R);
				G /= ctt;
				G = (G>255 ? 255 : G);
				G = (G<0 ? 0 : G);
				B /= ctt;
				B = (B>255 ? 255 : B);
				B = (B<0 ? 0 : B);

				o_data[o_pos++] = (unsigned char)B;
				o_data[o_pos++] = (unsigned char)G;
				o_data[o_pos++] = (unsigned char)R;

			}
		}
	}
	/*void DownSample(unsigned char* in, unsigned char* out, int width, int height, BayerType type)
	{
		int dwidth = width / 2;
		int dheight = height / 2;
		for (int i = 0; i < dheight; i++) {
			int line = i * dwidth;
			int offset = (i) % 2;
			for (int j = 0; j < dwidth; j++) {
				int ch = (j) % 2;
				int img_index = line + j;
				if (!offset) {
					switch (type)
					{
					case RGGB:
					case BGGR:
						if (!ch) {

							out[img_index] = in[i * 2 * width + j * 2];
						}
						else {
							out[img_index] = (in[i * 2 * width + j * 2 + 1] + in[(i * 2 + 1)* width + j * 2]) / 2; // G channel
						}
						break;
					case GRBG:
					case GBRG:
						if (!ch) {

							out[img_index] = (in[i * 2 * width + j * 2] + in[(i * 2 + 1)* width + j * 2 + 1]) / 2; // G channel
						}
						else {
							out[img_index] = in[i * 2 * width + j * 2 + 1];
						}
						break;
					default:
						//printf("error bayer type\n");
						return;
					}
				}
				else {

					switch (type)
					{
					case BGGR:
					case RGGB:
						if (!ch) {

							out[img_index] = (in[i * 2 * width + j * 2 + 1] + in[(i * 2 + 1)* width + j * 2]) / 2; // G channel
						}
						else {
							out[img_index] = in[(i * 2 + 1) * width + j * 2 + 1];
						}
						break;
					case GRBG:
					case GBRG:
						if (!ch) {
							out[img_index] = in[(i * 2 + 1)* width + j * 2];
						}
						else {
							out[img_index] = (in[i * 2 * width + j * 2] + in[(i * 2 + 1)* width + j * 2 + 1]) / 2; // G channel
						}
						break;
					default:
						//printf("error bayer type\n");
						return;
					}
				}

			}
		}
	}*/
}