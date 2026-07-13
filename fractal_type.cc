#include "fractal_type.hh"

#include <math.h>

inline void SixSide::per_image(double values[5], unsigned frame) const
{
}

inline void SixSide::per_line(double values[5], unsigned frame, unsigned line) const
{
}

inline void SixSide::per_pixel(double values[5], unsigned frame, unsigned line, unsigned row) const
{
}

std::vector<SixSide::SinCos> SixSide::calculate_sincos() const
{
	std::vector<SixSide::SinCos> result;
	result.reserve(c.size);
	for (unsigned i = 0; i < c.size; ++i)
	{
		const double phi = i*(3*M_PI)/c.size;
		const double ss = sin(phi);
		const double cc = cos(phi);
		result.emplace_back(fabs(ss),fabs(cc));
	}
	return result;
}