#ifndef IMAGERGBA_H
#define IMAGERGBA_H

#include <vector>
#include <string>
#include "Color.h"

struct RGBA{
	Byte r, g, b, a;
};

inline vec4 toVec4(RGBA C){
	return {
		ByteToFloat(C.r),
		ByteToFloat(C.g),
		ByteToFloat(C.b),
		ByteToFloat(C.a),
	};
}

class ImageRGBA{
	public:
	int width = 0, height = 0;
	std::vector<RGBA> pixels;

	ImageRGBA() = default;

	ImageRGBA(std::string filename);

	RGBA& operator()(int x, int y){
		return pixels[y*width + x];
	}
};

#endif
