
#ifndef FRACTAL_TYPE_HH
#define FRACTAL_TYPE_HH

#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Dimension
{
	unsigned size;
	double min;
	double max;
};

struct SinCos { double s, c; };


class FractalType
{
public:
	FractalType(const Dimension& _a, const Dimension& _b)
	: a(_a), b(_b)
	{}

	virtual ~FractalType() = default;

	virtual void per_image(double values[5], unsigned frame) const = 0;
	virtual void per_line (double values[5], unsigned frame, unsigned line) const = 0;
	virtual void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const = 0;
	[[nodiscard]] virtual const char* name() const = 0;

protected:
	const Dimension& a;
	const Dimension& b;
};


class SimpleFractalType : public FractalType
{
public:
	SimpleFractalType(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: FractalType(_a, _b), c(_c), d(_d), e(_e)
	{}

	void per_image(double values[5], unsigned frame) const override
	{
		values[2] = (c.size > 0) ? c.min + (c.max-c.min)/c.size*(frame+0.5) : c.min;
		values[3] = d;
		values[4] = e;
	}

	void per_line (double values[5], unsigned frame, unsigned line) const override
	{
		values[0] = a.min + (a.max-a.min)/a.size*(line+0.5);
	}

	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override
	{
		values[1] = b.min + (b.max-b.min)/b.size*(column+0.5);
	}

	[[nodiscard]] const char* name() const override { return "plain"; }

protected:
	const Dimension c;
	const double d;
	const double e;
};


class SixSide : public SimpleFractalType
{
	public:
	SixSide(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: SimpleFractalType(_a, _b, _c, _d, _e), sin_cos{calculate_sincos()}
	{}

	void per_image(double values[5], unsigned frame) const override;
	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "six_side"; }

private:

	[[nodiscard]] std::vector<SinCos> calculate_sincos() const;
	double sin_slide(double phi, const Dimension& dim) const
	{
		const double ds = dim.max - dim.min;
		return dim.max - ds*std::cos(phi);
	}

	const std::vector<SinCos> sin_cos;
};


class Circle : public SimpleFractalType
{
public:
	Circle(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	{}

	void per_image(double values[5], unsigned frame) const override;
	[[nodiscard]] const char* name() const override { return "circle"; }
};


class Sphere : public SimpleFractalType
{
public:
	Sphere(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, center_a((a.max+a.min)/2), radius_a((a.max-a.min)/2)
	, center_b((b.max+b.min)/2), radius_b((b.max-b.min)/2)
	, center_c((c.max+c.min)/2), radius_c((c.max-c.min)/2)
	, sc_line{ calculate_line_vector()}, sc_column{ calculate_column_vector()}
	, radius_de(_params.size()>1 ? _params[0] : 0.5)
	{}

	void per_image(double values[5], unsigned frame) const override;
	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;

private:
	[[nodiscard]] std::vector<SinCos> calculate_line_vector() const;
	[[nodiscard]] std::vector<SinCos> calculate_column_vector() const;

	const double center_a, radius_a, center_b, radius_b, center_c, radius_c;
	const std::vector<SinCos> sc_line;
	const std::vector<SinCos> sc_column;
	const double radius_de;
};


class Tileable : public SimpleFractalType
{
public:
	Tileable(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, center_a((a.max+a.min)/2), radius_a((a.max-a.min)/2)
	, center_b((b.max+b.min)/2), radius_b((b.max-b.min)/2)
	, sin_line{ calculate_sin_line()}, cos_column{ calculate_cos_column()}
	{}

	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "tileable"; }

private:
	[[nodiscard]] std::vector<double> calculate_sin_line() const;
	[[nodiscard]] std::vector<double> calculate_cos_column() const;

	const double center_a, radius_a, center_b, radius_b;
	const std::vector<double> sin_line;
	const std::vector<double> cos_column;
};


class TileableLoop : public Tileable
{
public:
	TileableLoop(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: Tileable(_a, _b, _c, _d, _e)
	{}

	void per_image(double values[5], unsigned frame) const override;
	[[nodiscard]] const char* name() const override { return "tileable-loop"; }
};

class Zoom : public SimpleFractalType
{
public:
	// _params: total zoom factor for the whole animation (final frame scale relative to frame 0),
	//          total rotation in degrees for the whole animation (default 0.1, 360)
	Zoom(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, center_a((a.max+a.min)/2), center_b((b.max+b.min)/2), center_c((c.max+c.min)/2)
	, zoom_total(_params.size()>0 ? _params[0] : 0.1)
	, rot_total(_params.size()>1 ? _params[1] * M_PI/180.0 : 2*M_PI)
	, cache_frame(~0u), a_x(0.0), a_y(0.0), b_x(0.0), b_y(0.0)
	{}

	void per_image(double values[5], unsigned frame) const override;
	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "zoom"; }

private:
	const double center_a, center_b, center_c;
	const double zoom_total, rot_total;
	// per-frame rotation/scale coefficients, cached by per_image (called once per frame before the parallel pixel loop)
	mutable unsigned cache_frame;
	mutable double a_x, a_y, b_x, b_y;
};
#endif // FRACTAL_TYPE_HH
