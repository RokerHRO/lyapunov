
#ifndef FRACTAL_TYPE_HH
#define FRACTAL_TYPE_HH
#include <string_view>

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

	virtual void per_image(double values[5]) = 0;
	virtual void per_line (double values[5], unsigned line) = 0;
	virtual void per_pixel(double values[5], unsigned line, unsigned row) = 0;

protected:
	const Dimension& a;
	const Dimension& b;
};


class SimpleFractalType : public FractalType
{
public:
	SimpleFractalType(const Dimension& _a, const Dimension& _b, double _c, double _d, double _e)
	: FractalType(_a, _b), c(_c), d(_d), e(_e)
	{}

	void per_image(double values[5]) override
	{
		values[2] = c;
		values[3] = d;
		values[4] = e;
	}

	void per_line (double values[5], unsigned line) override
	{
		values[0] = a.min + (a.max-a.min)/a.size*(line+0.5);
	}

	void per_pixel(double values[5], unsigned line, unsigned row) override
	{
		values[1] = b.min + (b.max-b.min)/b.size*(row+0.5);
	}

private:
	const double c;
	const double d;
	const double e;
};

#endif // FRACTAL_TYPE_HH