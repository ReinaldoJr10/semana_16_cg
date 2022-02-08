#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <string>
#include "Color.h"

class Image{
	public:
	int width = 0, height = 0;
	std::vector<Color> pixels;

	// Construtor padrão, imagem vazia
	Image() = default;

	// Construtor com tamanho da imagem
	Image(int w, int h):
		width(w), height(h), pixels(w*h)
	{}

	// Construtor para carregar um arquivo de imagem
	Image(std::string filename);

	// Acesso a um pixel
	Color& operator()(int x, int y){
		return pixels[y*width + x];
	}

	// Preenche a imagem com uma cor
	void fill(Color color){
		pixels.assign(pixels.size(), color);
	}

	// Salva a imagem num arquivo PNG
	void savePNG(std::string filename);

	// Salva a imagem num arquivo BMP
	void saveBMP(std::string filename);
};


#endif
