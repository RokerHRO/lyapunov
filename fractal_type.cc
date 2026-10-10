#include "fractal_type.hh"

#include <math.h>
#include <algorithm>

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
	const double phi = (c.size > 0) ? frame * 2*M_PI / c.size : 0.0;
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

void Tileable::per_line(double values[5], unsigned frame, unsigned line) const
{
	values[0] = center_a + radius_a * sin_line.at(line);
}

void Tileable::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	values[1] = center_b + radius_b * cos_column.at(column);
}

std::vector<double> Tileable::calculate_sin_line() const
{
	std::vector<double> result;
	result.reserve(a.size);
	for (unsigned i = 0; i < a.size; ++i)
		result.push_back(sin(i*(2*M_PI)/a.size));
	return result;
}

std::vector<double> Tileable::calculate_cos_column() const
{
	std::vector<double> result;
	result.reserve(b.size);
	for (unsigned i = 0; i < b.size; ++i)
		result.push_back(cos(i*(2*M_PI)/b.size));
	return result;
}

void TileableLoop::per_image(double values[5], unsigned frame) const
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

void Zoom::per_image(double values[5], unsigned frame) const
{
	if (frame != cache_frame)
	{
		const double span = (c.size > 1) ? (double)(c.size - 1) : 1.0;
		const double t = (double)frame / span;
		const double s = pow(zoom_total, t);
		const double theta = rot_total * t;
		a_x = (a.max-a.min) * s * sin(theta);
		a_y = (a.max-a.min) * s * cos(theta);
		b_x = (b.max-b.min) * s * cos(theta);
		b_y = (b.max-b.min) * s * sin(theta);
		cache_frame = frame;
	}
	values[2] = center_c;
	values[3] = d;
	values[4] = e;
}

void Zoom::per_line(double values[5], unsigned frame, unsigned line) const
{
}

void Zoom::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const double x = (column + 0.5) / b.size - 0.5;
	const double y = (line + 0.5) / a.size - 0.5;
	values[0] = center_a + a_x * x + a_y * y;
	values[1] = center_b + b_x * x - b_y * y;
}

void Wavy::per_line(double values[5], unsigned frame, unsigned line) const
{

}

void Wavy::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5)   + amp_a * sin(k_a * column - phase_a);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5) + amp_b * sin(k_b * line - phase_b);
}

detail::Perlin::Perlin(unsigned seed)
{
	int p[256];
	for (int i = 0; i < 256; ++i)
		p[i] = i;
	unsigned int s = seed ? seed : 1;
	for (int i = 255; i > 0; --i)
	{
		s = s*1103515245u + 12345u;
		const int j = (int)((s >> 16) % (unsigned)(i+1));
		std::swap(p[i], p[j]);
	}
	for (int i = 0; i < 512; ++i)
		perm[i] = p[i & 255];
}

double detail::Perlin::noise(double x, double y) const
{
	const int xi = (int)floor(x) & 255;
	const int yi = (int)floor(y) & 255;
	const double xf = x - floor(x);
	const double yf = y - floor(y);
	const double u = xf*xf*xf*(xf*(xf*6-15)+10);
	const double v = yf*yf*yf*(yf*(yf*6-15)+10);

	auto dot2 = [](int hash, double dx, double dy) -> double
	{
		switch (hash & 7)
		{
		case 0: return  dx + dy;
		case 1: return  dx - dy;
		case 2: return -dx + dy;
		case 3: return -dx - dy;
		case 4: return  dx;
		case 5: return -dx;
		case 6: return  dy;
		default: return -dy;
		}
	};

	const double aa = dot2(perm[perm[xi] + yi],     xf,     yf);
	const double ba = dot2(perm[perm[xi+1] + yi],   xf-1,   yf);
	const double ab = dot2(perm[perm[xi] + yi+1],   xf,     yf-1);
	const double bb = dot2(perm[perm[xi+1] + yi+1], xf-1,   yf-1);
	const double x1 = aa + u*(ba - aa);
	const double x2 = ab + u*(bb - ab);
	return x1 + v*(x2 - x1);
}

void Noise::per_line(double values[5], unsigned frame, unsigned line) const
{
	// A depends on both line and column (2D noise)
}

void Noise::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const double x = (column + 0.5) / b.size * freq;
	const double y = (line + 0.5) / a.size * freq;
	const double t = (c.size > 1) ? (double)frame / (c.size - 1) : 1.0;
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5) + amp_a * t * perlin.noise(x, y);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5) + amp_b * t * perlin.noise(x, y + 52.0);
}

void Inversion::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	// reflection about the inversion circle: distance d and angle theta from the center,
	// projected to the radius R^2/d (radius R scaled by R/d) at the same angle
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5);
	const double da = values[0] - center_a;
	const double db = values[1] - center_b;
	const double d = sqrt(da*da + db*db);
	if (d > 1e-150) // the center itself is fixed (its image is at infinity)
	{
		const double scale = radius*radius / (d*d);
		values[0] = center_a + scale * da;
		values[1] = center_b + scale * db;
	}
}
