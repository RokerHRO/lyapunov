
#ifndef LISSAJOUS_HH
#define LISSAJOUS_HH

#include <vector>
#include <ostream>

constexpr unsigned NumberOfPoints = 1024;

class Lissajous {
public:
	struct XY { double x,y; };

	Lissajous(unsigned _num, unsigned _denom, double _phase = 0.0)
	: num(_num), denom(_denom), phase(_phase)
	, points(generate_points())
	, length(generate_length())
	{}

	void gen_svg(std::ostream& os) const;

	const unsigned num;
	const unsigned denom;
	const double phase;
	const std::vector<XY> points;
	const double length;

	[[nodiscard]] std::vector<XY> generate_points() const;
	[[nodiscard]] double generate_length() const;
};

#endif // LISSAJOUS_HH
