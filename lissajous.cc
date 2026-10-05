
#include "lissajous.hh"

#include <cmath>

void Lissajous::gen_svg(std::ostream& os) const {
	XY p0 = points.at(0);
	os << "<g id=\"lissa_" << num << "_" << denom << "_" << phase << "\" >\n"
		"\t<rect x=\"0\" y=\"0\" width=\"256\" height=\"256\" fill=\"#ffeedd\"/>\n"
		"\t<path fill=\"none\" stroke=\"blue\" d=\"M" << p0.x*128+128 << " " << p0.y*128+128;
	for (unsigned u=1; u<points.size(); ++u)
	{
		const XY p = points.at(u);
		os << " L" << p.x*128 + 128 << " " << p.y*128+128;
	}
	os << "\" />\n</g>";
}

std::vector<Lissajous::XY> Lissajous::generate_points() const {
	std::vector<XY> vp;
	vp.reserve(NumberOfPoints);
	for (unsigned int i = 0; i < NumberOfPoints; ++i)
	{
		const double phi = i * (2*M_PI/NumberOfPoints);
		vp.emplace_back(XY{sin(phi*num), cos(phi*denom + phase)} );
	}
	return vp;
}

double Lissajous::generate_length() const {
	if (points.empty())
		return 0.0;

	double len = 0;
	XY p0 = points.at(0);
	for (int i = 1; i < NumberOfPoints; ++i)
	{
		const XY p1 = points.at(i);
		len += hypot(p1.x - p0.x, p1.y - p0.y);
		p0 = p1;
	}
	return len;
}
