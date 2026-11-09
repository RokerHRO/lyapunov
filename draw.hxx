#ifndef DRAW_HXX
#define DRAW_HXX

#include "draw.hh"

#include <algorithm>
#include <cstdlib>

template <class PXL>
void draw_hline(std::vector<PXL>& img, const int width, const int height, const int x0, const int x1, const int y, const PXL& color)
{
	if (y < 0 || y >= height)
		return;
	const int xa = std::max(std::min(x0, x1), 0);
	const int xb = std::min(std::max(x0, x1), width-1);
	if (xa > xb)
		return;
	PXL* p = &img[size_t(y) * width + xa];
	for (int i = 0; i <= xb - xa; ++i)
		p[i] = color;
}

template <class PXL>
void draw_vline(std::vector<PXL>& img, const int width, const int height, const int x, const int y0, const int y1, const PXL& color)
{
	if (x < 0 || x >= width)
		return;
	const int ya = std::max(std::min(y0, y1), 0);
	const int yb = std::min(std::max(y0, y1), height-1);
	if (ya > yb)
		return;
	PXL* p = &img[size_t(ya) * width + x];
	for (int i = 0; i <= yb - ya; ++i)
		p[i*width] = color;
}

template <class PXL>
void draw_line(std::vector<PXL>& img, const int width, const int height, int x0, int y0, const int x1, const int y1, const PXL& color)
{
	if (y0 == y1)
	{
		draw_hline(img, width, height, x0, x1, y0, color);
		return;
	}
	if (x0 == x1)
	{
		draw_vline(img, width, height, x0, y0, y1, color);
		return;
	}
	// integer Bresenham
	const int dx = std::abs(x1 - x0);
	const int dy = std::abs(y1 - y0);
	const int sx = x0 < x1 ? 1 : -1;
	const int sy = y0 < y1 ? 1 : -1;
	int err = dx - dy;
	for (;;)
	{
		if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height)
			img[size_t(y0) * width + x0] = color;
		if (x0 == x1 && y0 == y1)
			break;
		const int e2 = 2*err;
		if (e2 > -dy)
		{
			err -= dy;
			x0 += sx;
		}
		if (e2 < dx)
		{
			err += dx;
			y0 += sy;
		}
	}
}

template <class PXL>
void draw_text(std::vector<PXL>& /*img*/, const int /*width*/, const int /*height*/, const int /*x*/, const int /*y*/, const std::string& /*text*/, const PXL& /*fg*/, const PXL& /*bg*/, const int /*font_size*/)
{
}

#endif // DRAW_HXX
