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
void draw_character(std::vector<PXL>& img, const int width, const int height,
                    const int x, const int y,
                    const char c, const PXL& color, const Font& font)
{
	// data[] only covers 0x20..0x7F; characters outside the range fall back to 0x7F
	const unsigned char uc = static_cast<unsigned char>(c);
	const int index = uc < 0x20 || uc > 0x7f ? 0x7f - 0x20 : uc - 0x20;
	const uint16_t* glyph = &font.data[size_t(index) * font.height];
	for (unsigned ry = 0; ry < font.height; ++ry)
	{
		const int iy = y + int(ry);
		if (iy < 0 || iy >= height)
			continue;
		const uint16_t row = glyph[ry];
		// row is MSB-first: pixel rx (0 = leftmost) is bit font.width-1-rx
		for (unsigned rx = 0; rx < font.width; ++rx)
		{
			if ((row & (uint16_t(1) << (font.width - 1 - rx))) == 0)
				continue;
			const int ix = x + int(rx);
			if (ix < 0 || ix >= width)
				continue;
			img[size_t(iy) * width + ix] = color;
		}
	}
}


template <class PXL>
void draw_text(std::vector<PXL>& img, const int width, const int height,
               const int x, const int y,
               const std::string& text, const PXL& color, const Font& font)
{
	// '\n' continues at the original x, font.height further down
	const int x0 = x;
	int cx = x;
	int cy = y;
	for (const char c : text)
	{
		if (c == '\n')
		{
			cx = x0;
			cy += int(font.height);
			continue;
		}
		draw_character(img, width, height, cx, cy, c, color, font);
		cx += int(font.width);
	}
}

template <class PXL>
void draw_shaded_text(std::vector<PXL>& img, int width, int height, int x, int y, const std::string& text, const PXL& foreground, const PXL& shadow, const Font& font)
{
	draw_text(img, width, height, x+1, y, text, shadow, font);
	draw_text(img, width, height, x+1, y+1, text, shadow, font);
	draw_text(img, width, height, x, y+1, text, shadow, font);
	draw_text(img, width, height, x, y, text, foreground, font);
}

#endif // DRAW_HXX
