#ifndef DRAW_HH
#define DRAW_HH

#include <string>
#include <vector>

#include "rgb.hh"

// public interface: draw into a row-major pixel buffer of size_t(width) * height
// pixels, x = column (to the right), y = row (downwards)
template <class PXL>
void draw_hline(std::vector<PXL>& img, const int width, const int height, const int x0, const int x1, const int y, const PXL& color);

template <class PXL>
void draw_vline(std::vector<PXL>& img, const int width, const int height, const int x, const int y0, const int y1, const PXL& color);

template <class PXL>
void draw_line(std::vector<PXL>& img, const int width, const int height, int x0, int y0, const int x1, const int y1, const PXL& color);

template <class PXL>
void draw_text(std::vector<PXL>& img, const int width, const int height, const int x, const int y, const std::string& text, const PXL& fg, const PXL& bg, const int font_size);

#endif // DRAW_HH
