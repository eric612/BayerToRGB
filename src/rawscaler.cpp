#include "rawscaler.h"
#include "math.h"
#include "util.h"
#include "stdlib.h"
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef abs_min
#define abs_min( a, b ) ( ((abs(a)) < (abs(b))) ? (a) : (b) )
#endif
#ifndef abs_max
#define abs_max( a, b ) ( ((abs(a)) > (abs(b))) ? (a) : (b) )
#endif
#define BOUND(a,min_val,max_val)           ( (a < min_val) ? min_val : (a >= max_val) ? (max_val) : a )
#include <string.h>

#ifdef openmp
#include <omp.h>
#endif
namespace RawScaler {


	float RawProcessor::AdamsInterpolationOri(unsigned char* in, int x, int y, int width, float high_pass_weight = 1.0)
	{
		int v1 = in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x];
		int v2 = in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2];
		int green = in[(y + 1) * width + x] + in[(y - 1) * width + x] + in[y * width + x + 1] + in[y * width + x - 1];
		return (green + v1 + v2)/4.0;
	}
	float RawProcessor::AdamsInterpolation(unsigned char* in, int x, int y, int width,float high_pass_weight = 1.0)
	{
		float max_cut = 30;
		int v1 = in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x];
		int v2 = in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2];
		int direction = fabs(v2) < fabs(v1);
		if (direction == 0)
			//return  BOUND(((in[(y + 1) * width + x] + in[(y - 1) * width + x])*0.5 + (v1)*0.25), 1, 255);
			return  ((in[(y + 1) * width + x] + in[(y - 1) * width + x])*0.5 + BOUND((v1), -max_cut, max_cut)*high_pass_weight*0.25);
		else
			//return  BOUND(((in[y * width + x + 1] + in[y * width + x - 1])*0.5 + (v2)*0.25), 1, 255);
			return  ((in[y * width + x + 1] + in[y * width + x - 1])*0.5 + BOUND((v2), -max_cut, max_cut)*high_pass_weight*0.25);
	}

	float RawProcessor::AdamsInterpolation (unsigned char* in, int x, int y, int width, int direction, float max_cut = 30)
	{
		if(direction == 0)
			//return BOUND(((in[(y + 1) * width + x] + in[(y - 1) * width + x])*0.5 + (in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x])*0.25),1,255);
			return ((in[(y + 1) * width + x] + in[(y - 1) * width + x])*0.5 + BOUND((in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x]),-max_cut, max_cut)*0.25);
		else 
			//return BOUND(((in[y * width + x + 1] + in[y * width + x - 1])*0.5 + (in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2])*0.25), 1, 255);
			return ((in[y * width + x + 1] + in[y * width + x - 1])*0.5 + BOUND((in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2]), -max_cut, max_cut)*0.25);

	}

	float RawProcessor::Interpolation_G_Only(unsigned char* in, int x, int y, int width, float max_cut = 64,int kernel_size = 5)
	{
		float g_arr[8];
		float r, g, b;
		float mean[3][2];
		float *gbr = new float[3];
		int b1 = in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x];
		int b2 = in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2];
		int b3 = in[(y + 1) * width + x] + in[(y - 1) * width + x];
		int b4 = in[y * width + x + 1] + in[y * width + x - 1];


		
		float g0h = AdamsInterpolation(in, x, y, width, 1, max_cut);
		float g1h = AdamsInterpolation(in, x + 2, y, width, 1, max_cut);
		float g2h = AdamsInterpolation(in, x - 2, y, width, 1, max_cut);
		float g3h = AdamsInterpolation(in, x - 1, y - 1, width, 1, max_cut);
		float g4h = AdamsInterpolation(in, x - 1, y + 1, width, 1, max_cut);
		float g5h = AdamsInterpolation(in, x + 1, y - 1, width, 1, max_cut);
		float g6h = AdamsInterpolation(in, x + 1, y + 1, width, 1, max_cut);

		float g0v = AdamsInterpolation(in, x, y, width, 0, max_cut);
		float g1v = AdamsInterpolation(in, x, y + 2, width, 0, max_cut);
		float g2v = AdamsInterpolation(in, x, y - 2, width, 0, max_cut);
		float g3v = AdamsInterpolation(in, x - 1, y - 1, width, 0, max_cut);
		float g4v = AdamsInterpolation(in, x - 1, y + 1, width, 0, max_cut);
		float g5v = AdamsInterpolation(in, x + 1, y - 1, width, 0, max_cut);
		float g6v = AdamsInterpolation(in, x + 1, y + 1, width, 0, max_cut);

		float r0h = AdamsInterpolation(in, x + 1, y, width, 1, max_cut);
		float r1h = AdamsInterpolation(in, x - 1, y, width, 1, max_cut);
		float r0v = AdamsInterpolation(in, x, y + 1, width, 0, max_cut);
		float r1v = AdamsInterpolation(in, x, y - 1, width, 0, max_cut);
		r = in[y * width + x];
		float grh[9], grv[9], rh[9], rv[9];
		float maxima_diff = 1.0;
		grh[0] = (g0h - in[y * width + x]);
		grh[1] = (in[y * width + x + 1] - r0h);

		grh[2] = (in[y * width + x - 1] - r1h);
		grh[3] = (g1h - in[y * width + x + 2]);
		grh[4] = (g2h - in[y * width + x - 2]);
		grh[5] = g3h;
		grh[6] = g4h;
		grh[7] = g5h;
		grh[8] = g6h;

		grv[0] = (g0v - in[y * width + x]);
		grv[1] = (in[(y + 1) * width + x] - r0v);

		grv[2] = (in[(y - 1) * width + x] - r1v);
		grv[3] = (g1v - in[(y + 2) * width + x]);
		grv[4] = (g2v - in[(y - 2) * width + x]);
		grv[5] = g3v;
		grv[6] = g4v;
		grv[7] = g5v;
		grv[8] = g6v;

		rh[0] = in[y * width + x];
		rh[1] = r0h;
		rh[2] = r1h;
		rh[3] = in[y * width + x + 2];
		rh[4] = in[y * width + x - 2];

		rv[0] = in[y * width + x];
		rv[1] = r0v;
		rv[2] = r1v;
		rv[3] = in[(y + 2) * width + x];
		rv[4] = in[(y - 2) * width + x];

		float gr_std_h = GetStd(grh, kernel_size);
		float gr_std_v = GetStd(grv, kernel_size);

		float alpha1 = 0.5;
		float alpha2 = 0.5;
		if (gr_std_h + gr_std_v != 0) {
			float tmp1 = pow(gr_std_h, 3.0);
			float tmp2 = pow(gr_std_v, 3.0);
			if (tmp1 + tmp2 > 0) {
				alpha1 = tmp1 / (tmp1 + tmp2);
				alpha2 = tmp2 / (tmp1 + tmp2);
			}

		}

		g0v = (in[(y + 1) * width + x] + in[(y - 1) * width + x])*0.5 + (in[y * width + x] * 2 - in[(y + 2) * width + x] - in[(y - 2) * width + x])*0.25 * 0;
		g0h = (in[y * width + x + 1] + in[y * width + x - 1])*0.5 + (in[y * width + x] * 2 - in[y  * width + x + 2] - in[y  * width + x - 2])*0.25 * 0;

		float certain_h = (in[y * width + x] * 1) / 1;
		float certain_v = (in[y * width + x] * 1) / 1;
		float tmp2 = (b3 * 2 + BOUND(b1 * 1, -10, 10)) / 4;
		float tmp1 = (b4 * 2 + BOUND(b2 * 1, -10, 10)) / 4;
		tmp1 = BOUND((certain_h + GetMedian(grh, 3)),0,255);
		tmp2 = BOUND((certain_v + GetMedian(grv, 3)),0,255);
		g_arr[0] = tmp1*alpha2 + tmp2*alpha1;

		g = g_arr[0];
		if (g > 255) g = 255;
		else if (g < 0) g = 0;
		delete[] gbr;
		return g;
	}


	int RawProcessor::GetCurrentColorType(int offset_y,int offset_x,BayerType type)
	{
		enum {
			B = 0,
			G = 1,
			R = 2
		};
		switch (type)
		{
		case RawScaler::RGGB:
			break;
		case RawScaler::GRBG:
			break;
		case RawScaler::BGGR:
			if (!offset_y) {
				if (!offset_x) {
					return B;
				}
				else {
					return G;
				}
			}
			else {
				if (!offset_x) {
					return G;
				}
				else {
					return R;
				}
			}
			break;
		case RawScaler::GBRG:
			break;
		default:
			break;
		}
	}
	void Scale_Down_By_Area_only_Y(int o_width, int o_height, unsigned char *o_data, int i_width, int i_height, unsigned char *i_data)
	{
		int st_x, end_x;
		int st_y, end_y;
		int x, y;
		int y1, x1;
		int i_pos;
		int o_pos;
		int w_line_pixel;
		double ctt;
		double weight;

		int Y;
		double w1, h1;
		double x1_pos;
		double y1_pos;
		double x1_left;
		double x1_right;
		double y1_top;
		double y1_bottom;
		double radius_w;
		double radius_h;

		double dx, dy;
		double dx_left;
		double dx_right;
		double dy_top;
		double dy_bottom;

		//input image one pixel's w=1 and h=1
		//output image one pixel's w1 and h1
		//w*i_width = w1*o_width;
		w1 = (double)i_width / (double)o_width;
		h1 = (double)i_height / (double)o_height;

		radius_w = w1 / 2;
		radius_h = h1 / 2;

		w_line_pixel = i_width;
		o_pos = 0;
		ctt = w1*h1;

		for (y1 = 0; y1<o_height; y1++)
		{

			y1_pos = (y1 + 0.5)*h1;
			y1_top = y1_pos - radius_h;
			y1_bottom = y1_pos + radius_h;

			st_y = int(y1_top);
			end_y = int(y1_bottom + 0.9999) - 1;

			dy_top = 1 - (y1_top - st_y);
			dy_bottom = y1_bottom - end_y;

			for (x1 = 0; x1<o_width; x1++)
			{
				x1_pos = (x1 + 0.5)*w1;
				//center (x1_pos,y1_pos)
				//calculate area point
				x1_left = x1_pos - radius_w;
				x1_right = x1_pos + radius_w;

				st_x = int(x1_left);
				end_x = int(x1_right + 0.9999) - 1;

				dx_left = 1 - (x1_left - st_x);
				dx_right = x1_right - end_x;

				Y = 0;
				//ctt=0;

				//calculate area
				for (y = st_y; y <= end_y; y++)
				{
					if (y == st_y)
						dy = dy_top;
					else if (y == end_y)
						dy = dy_bottom;
					else
						dy = 1;

					i_pos = y*w_line_pixel + st_x;

					for (x = st_x; x <= end_x; x++)
					{
						if (x == st_x)
							dx = dx_left;
						else if (x == end_x)
							dx = dx_right;
						else
							dx = 1;

						weight = dx*dy;
						Y += (int)(weight*i_data[i_pos++]);
						//ctt+=weight;
					}
				}

				Y = (int)(Y / ctt);
				Y = (Y>255 ? 255 : Y);

				o_data[o_pos++] = (unsigned char)Y;


			}
		}
	}
	void Scale_Up_By_Bilinear_only_Y(int o_width, int o_height, unsigned char *o_data, int i_width, int i_height, unsigned char *i_data)
	{

		int x, y;
		int y1, x1;
		double w1, h1;
		double x1_pos;
		double y1_pos;
		int dx, dy;
		double x1_offset;
		double y1_offset;


		int *x_hash;
		int *dx_hash;

		x_hash = new int[o_width];
		dx_hash = new int[o_width];


		//input image one pixel's w=1 and h=1
		//output image one pixel's w1 and h1
		//w*i_width = w1*o_width;
		w1 = (double)i_width / (double)o_width;
		h1 = (double)i_height / (double)o_height;

		x1_offset = w1 / 2.0;
		y1_offset = h1 / 2.0;

		int Y;
		int weight;
		int o_y_pos = 0;
		int i_y_pos = 0;


		for (x1 = 0; x1<o_width; x1++)
		{

			x1_pos = x1*w1 + x1_offset;
			x = int(x1_pos - 0.5);
			dx = (int)(128 * (x1_pos - x - 0.5));
			if (dx <= 0)
				dx = 0;
			if (x >= i_width)
				x = i_width - 1;

			x_hash[x1] = x;
			dx_hash[x1] = dx;
		}

		int p[4];
		o_y_pos = 0;
		for (y1 = 0; y1<o_height; y1++)
		{

			y1_pos = y1*h1 + y1_offset;
			y = int(y1_pos - 0.5);
			dy = (int)(128 * (y1_pos - y - 0.5));
			if (y >= i_height)
				y = i_height - 1;
			if (dy <= 0)
				dy = 0;
			for (x1 = 0; x1<o_width; x1++)
			{


				x = x_hash[x1];
				dx = dx_hash[x1];
				i_y_pos = y*i_width + x;
				Y = 0;

				p[0] = i_data[i_y_pos];
				if (x<i_width - 1)
					p[1] = i_data[i_y_pos + 1];
				else
					p[1] = p[0];

				if (y<i_height - 1)
				{
					p[2] = i_data[i_y_pos + i_width];
					if (x<i_width - 1)
						p[3] = i_data[i_y_pos + i_width + 1];
					else
						p[3] = p[2];
				}
				else
				{
					p[2] = p[0];
					p[3] = p[1];
				}

				weight = (128 - dx)*(128 - dy);
				Y = p[0] * weight;
				weight = (dx)*(128 - dy);
				Y += p[1] * weight;
				weight = (128 - dx)*(dy);
				Y += p[2] * weight;
				weight = (dx)*(dy);
				Y = (Y + p[3] * weight);
				Y = Y >> 14;
				if (Y>255)
					Y = 255;
				if (Y<0)
					Y = 0;

				o_data[o_y_pos++] = (unsigned char)Y;

			}


		}
		delete[] x_hash;
		delete[] dx_hash;

	}
	int RawProcessor::GetValFromArray(int x,int y,unsigned char* data,int width)
	{
		return data[y*width + x];
	}
	void RawProcessor::BayerToRGB(unsigned char* in, unsigned char* out, int iwidth, int iheight, int owidth, int oheight, BayerType type)
	{
		float w_ratio = iwidth / (float)owidth;
		float h_ratio = iheight / (float)oheight;

		int dwidth = owidth;
		int dheight = oheight;
		float ratio = w_ratio * h_ratio;
		int rb_width = owidth;
		int rb_height = oheight;
		if (w_ratio < 2 || h_ratio < 2) {
			//rb_width /= 2;
			//rb_height /= 2;
		}
		unsigned char* temp_g = new unsigned char[iwidth*iheight];
		unsigned char* temp_gr = new unsigned char[iwidth*iheight];
		unsigned char* temp_gb = new unsigned char[iwidth*iheight];
		unsigned char* temp_g2 = new unsigned char[owidth*oheight];
		unsigned char* temp_g2r = new unsigned char[rb_width*rb_height];
		unsigned char* temp_g2b = new unsigned char[rb_width*rb_height];
#ifdef openmp
#pragma omp parallel for
#endif
		for (int i = 4; i < iheight - 4; i++) {
			int line = i * iwidth;
			int offset = (i) % 2;
			//dy = dy - (dy % 2);
			for (int j = 0; j < iwidth - 0; j++) {
				int ch = (j) % 2;
				int img_index = line + j;
				int img_rb_index = (i/2) * (iwidth/2) + j/2;
				int val[3];
				switch (type)
				{
				case RGGB:
					if (!offset) {
						if (ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gr[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					else {
						if (!ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gb[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					break;
				case BGGR:
					if (!offset) {
						if (ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							if (1) {
								temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);								
								temp_gb[img_rb_index] = BOUND(in[img_index] - temp_g[img_index],-255, 255) / 2 +127;
							}
						}
					}
					else {
						if (!ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							if (1) {
	
								temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
								temp_gr[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
							}
						}
					}

					break;
				case GRBG:
					if (!offset) {
						if (!ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gr[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					else {
						if (ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gb[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					break;
				case GBRG:
					if (!offset) {
						if (!ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gb[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					else {
						if (ch) {
							temp_g[img_index] = in[img_index];
						}
						else {
							temp_g[img_index] = Interpolation_G_Only(in, j, i, iwidth);
							temp_gr[img_rb_index] = BOUND(in[img_index] - temp_g[img_index], -255, 255) / 2 + 127;
						}
					}
					break;
				default:
					temp_g[img_index] = in[img_index];
					//return raw_img;
				}

			}
		}
		Scale_Down_By_Area_only_Y(owidth, oheight, temp_g2, iwidth, iheight, temp_g);
		if (w_ratio < 2 || h_ratio < 2) {
			Scale_Up_By_Bilinear_only_Y(rb_width, rb_height, temp_g2r, iwidth/2 , iheight / 2, temp_gr);
			Scale_Up_By_Bilinear_only_Y(rb_width, rb_height, temp_g2b, iwidth / 2, iheight / 2, temp_gb);
		}
		else {
			Scale_Down_By_Area_only_Y(rb_width, rb_height, temp_g2r, iwidth/2 , iheight/2 , temp_gr);
			Scale_Down_By_Area_only_Y(rb_width, rb_height, temp_g2b, iwidth/2 , iheight/2 , temp_gb);
		}

#ifdef openmp
#pragma omp parallel for
#endif
		for (int i = 2; i < dheight - 2; i++) {
			int line = i * dwidth;
			int offset = (i) % 2;
			int dy = floor(i*h_ratio);
			float fy = i*h_ratio - dy;
			//dy = dy - (dy % 2);
			for (int j = 0; j < dwidth - 0; j++) {
				int ch = (j) % 2;
				int img_index = line + j;
				int img_rb_index = (i / 2) * (dwidth / 2) + j / 2;
				if (w_ratio < 2 || h_ratio < 2) {
					img_rb_index = (i ) * (dwidth ) + j ;
				}

				
				//dy = dy - (dy % 2);
				int dx = floor(j*w_ratio);
				float fx = j*w_ratio - dx;
				int val[3];

				int bgr = GetCurrentColorType(offset, ch, type);
				int rb[2];
				out[img_index * 3 + 1] = temp_g2[img_index];
				if (1)
				{
					if (w_ratio < 2 || h_ratio < 2) {
						float bg[9];
						bg[0] = temp_g2b[img_rb_index - 1];
						bg[1] = temp_g2b[img_rb_index];
						bg[2] = temp_g2b[img_rb_index + 1];
						bg[3] = temp_g2b[img_rb_index - (dwidth / 2) - 1];
						bg[4] = temp_g2b[img_rb_index - (dwidth / 2)];
						bg[5] = temp_g2b[img_rb_index - (dwidth / 2) + 1];
						bg[6] = temp_g2b[img_rb_index + (dwidth / 2) - 1];
						bg[7] = temp_g2b[img_rb_index + (dwidth / 2)];
						bg[8] = temp_g2b[img_rb_index + (dwidth / 2) + 1];

						float rg[9];
						rg[0] = temp_g2r[img_rb_index - 1];
						rg[1] = temp_g2r[img_rb_index];
						rg[2] = temp_g2r[img_rb_index + 1];
						rg[3] = temp_g2r[img_rb_index - (dwidth / 2) - 1];
						rg[4] = temp_g2r[img_rb_index - (dwidth / 2)];
						rg[5] = temp_g2r[img_rb_index - (dwidth / 2) + 1];
						rg[6] = temp_g2r[img_rb_index + (dwidth / 2) - 1];
						rg[7] = temp_g2r[img_rb_index + (dwidth / 2)];
						rg[8] = temp_g2r[img_rb_index + (dwidth / 2) + 1];

						out[img_index * 3] = BOUND(temp_g2[img_index] + (temp_g2b[img_rb_index] - 127)*2, 0, 255);
						//out[img_index * 3] = BOUND(temp_g2[img_index] + (GetMedian(bg, 2) - 127), 0, 255);
						out[img_index * 3 + 2] = BOUND(temp_g2[img_index] + (temp_g2r[img_rb_index] - 127) * 2, 0, 255);
						//out[img_index * 3 + 2] = BOUND(temp_g2[img_index] + (GetMedian(rg, 2) - 127), 0, 255);
					}
					else {
						out[img_index * 3] = BOUND(temp_g2[img_index] + (temp_g2b[img_index] - 127) * 2, 0, 255);
						out[img_index * 3 + 2] = BOUND(temp_g2[img_index] + (temp_g2r[img_index] - 127) * 2, 0, 255);
					}

					
				}			
			}
		}
		delete[] temp_g;
		delete[] temp_g2;
		delete[] temp_gr;
		delete[] temp_gb;
		delete[] temp_g2r;
		delete[] temp_g2b;
	}

}
