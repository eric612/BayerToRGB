#ifndef RAWSCALER_H_
#define RAWSCALER_H_
#include "util.h"
namespace RawScaler {

	class RawProcessor {
	public:
		void BayerToRGB(unsigned char* in, unsigned char* out, int iwidth, int iheight, int owidth, int oheight, BayerType type);
	private:
		int GetValFromArray(int x, int y, unsigned char* data, int width);
		int GetCurrentColorType(int offset_y, int offset_x, BayerType type);
		float AdamsInterpolationOri(unsigned char* in, int x, int y, int width, float high_pass_weight );
		float AdamsInterpolation(unsigned char* in, int x, int y, int width, float high_pass_weight );
		float AdamsInterpolation(unsigned char* in, int x, int y, int width, int direction, float max_cut );
		float Interpolation_G_Only(unsigned char* in, int x, int y, int width, float max_cut, int kernel_size);
	};
}
#endif  // RAWSCALER_H_