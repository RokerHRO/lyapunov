#ifndef DRAW_HH
#define DRAW_HH

#include <cstdint>
#include <map>
#include <string>
#include <vector>


struct Font
{
	unsigned width; // in pixel. max. 16 pixel supported
	unsigned height;
	const char* name;
	const uint16_t* data; // contains only data for character codes 0x20 ... 0x7F.  'height' entries per character cell
};

// defined in fonts/font.cc
extern const std::map<std::string, const Font*> fonts;

// public interface: draw into a row-major pixel buffer of size_t(width) * height
// pixels, x = column (to the right), y = row (downwards)
template <class PXL>
void draw_hline(std::vector<PXL>& img, int width, int height, int x0, int x1, int y, const PXL& color);

template <class PXL>
void draw_vline(std::vector<PXL>& img, int width, int height, int x, int y0, int y1, const PXL& color);

template <class PXL>
void draw_line(std::vector<PXL>& img, int width, int height, int x0, int y0, int x1, int y1, const PXL& color);

template <class PXL>
void draw_text(std::vector<PXL>& img, int width, int height, int x, int y, const std::string& text, const PXL& color, const Font& font);

#endif // DRAW_HH
