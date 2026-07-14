
#ifndef FRACTAL_TYPE_HH
#define FRACTAL_TYPE_HH

#include <vector>

struct Dimension
{
	unsigned size;
	double min;
	double max;
};

class FractalType
{
public:
	FractalType(const Dimension& _a, const Dimension& _b)
	: a(_a), b(_b)
	{}

	virtual ~FractalType() = default;

	virtual void per_image(double values[5], unsigned frame) const = 0;
	virtual void per_line (double values[5], unsigned frame, unsigned line) const = 0;
	virtual void per_pixel(double values[5], unsigned frame, unsigned line, unsigned row) const = 0;
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
		values[2] = c.min + (c.max-c.min)/c.size*(frame+0.5);
		values[3] = d;
		values[4] = e;
	}

	void per_line (double values[5], unsigned frame, unsigned line) const override
	{
		values[0] = a.min + (a.max-a.min)/a.size*(line+0.5);
	}

	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned row) const override
	{
		values[1] = b.min + (b.max-b.min)/b.size*(row+0.5);
	}

	const char* name() const override { return "plain"; }

protected:
	const Dimension c;
	const double d;
	const double e;
};

class SixSide : public SimpleFractalType {
	public:
	SixSide(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: SimpleFractalType(_a, _b, _c, _d, _e), sin_cos{calculate_sincos()}
	{}

	void per_image(double values[5], unsigned frame) const override;
	void per_line(double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned row) const override;
	const char* name() const override { return "six_side"; }

private:
	struct SinCos {double s,c;};

	std::vector<SinCos> calculate_sincos() const;

	const std::vector<SinCos> sin_cos;
};


class Circle : public SimpleFractalType {
public:
	Circle(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	{}

	void per_image(double values[5], unsigned frame) const override;
	const char* name() const override { return "circle"; }

};

#endif // FRACTAL_TYPE_HH
