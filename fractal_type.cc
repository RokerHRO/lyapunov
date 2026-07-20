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

std::vector<SinCos> SixSide::calculate_sincos() const
{
	std::vector<SinCos> result;
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

void Circle::per_image(double values[5], unsigned frame) const
{
	if (c.size==0)
	{
		values[2] = c.min;
		values[3] = d;
		values[4] = e;
	}else{
		const double cm = (c.min + c.max)/2;
		const double r = abs(c.max - c.min);
		const double phi = frame * 2*M_PI / c.size;
		values[2] = cm - r*cos(phi);
		values[3] = d + r*sin(phi);
		values[4] = e;
	}
}


void Sphere::per_image(double values[5], unsigned frame) const
{
	const double phi = frame * 2*M_PI / c.size;
	values[3] = d - radius_de * cos(phi);
	values[4] = e + radius_de * sin(phi);
}

void Sphere::per_line(double values[5], unsigned frame, unsigned line) const
{
	values[0] = center_a - sc_line.at(line).c * radius_a;
}

void Sphere::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const auto pl = sc_line.at(line);
	const auto pc = sc_column.at(column);
	values[1] = center_b + pl.s * pc.c * radius_b;
	values[2] = center_b - pl.s * pc.s * radius_c;
}

std::vector<SinCos> Sphere::calculate_line_vector() const
{
	std::vector<SinCos> result;
	result.reserve(a.size);
	for (unsigned i = 0; i < a.size; ++i)
	{
		const double phi = i*M_PI/a.size;
		const double ss = sin(phi);
		const double cc = cos(phi);
		result.emplace_back(ss,cc);
	}
	return result;
}

std::vector<SinCos> Sphere::calculate_column_vector() const
{
	std::vector<SinCos> result;
	result.reserve(b.size);
	for (unsigned i = 0; i < b.size; ++i)
	{
		const double phi = i*(2*M_PI)/b.size;
		const double ss = sin(phi);
		const double cc = cos(phi);
		result.emplace_back(ss,cc);
	}
	return result;
}
