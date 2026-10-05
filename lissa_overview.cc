#include "lissajous.hh"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

constexpr unsigned XSpace = 15;
constexpr unsigned YSpace = 60;

struct ND
{
	unsigned n;
	unsigned d;
	unsigned ratio() const { return unsigned( n*2520.0/d + 0.5 ); }
};

bool operator<(const ND lhs, const ND rhs) { return lhs.ratio() < rhs.ratio(); }
bool operator==(const ND lhs, const ND rhs) { return lhs.ratio() == rhs.ratio(); }

std::vector<ND> nd;

int main()
{
	for (unsigned n=1; n<=8; ++n)
		for (unsigned d=1; d<=n; ++d)
		{
			ND value{d,n};
			if (std::ranges::find(nd, value) == nd.end())
			{
				nd.push_back(value);
			}
		}
	const unsigned Width = 7*(256+XSpace) + 66;
	const unsigned Height = (256+YSpace)*nd.size();
	std::cout << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
		"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << Width << "\" height=\"" << Height << "\" viewBox=\"0 0 " << Width << " " << Height << "\" >\n";
	std::cout << "<style> text { font-family:Arial, Helvetica, sans-serif; font-size:25px; fill:green; stroke:none;} </style>\n";
	std::cout << R"(<g stroke="black" stroke-width="2px">)";
	for (unsigned row=0; row<nd.size(); ++row)
	{
		const auto [n, d] = nd.at(row);
		std::cout << "<text x=\"0\" y=\"" << (row+0.5)*(256+YSpace) << "\">" << n << " : " << d << "</text>\n";
		for (int column=0; column <= 6; ++column)
		{
			Lissajous lissi(n, d, column*15*M_PI/180);
			std::cerr << "#" << row << "/" << column << ": n=" << n << ", d=" << d << " len=" << lissi.length << ". " << std::endl;
			std::cout << "<g transform=\"translate(" << column*(256+XSpace)+66 << ", " << row*(256+YSpace) << ")\" >\n";
			lissi.gen_svg(std::cout);
			std::cout << "<text x=\"0\" y=\"285\">length≈" << lissi.length << "</text>\n";
			std::cout << "</g>\n";
		}
	}
	std::cout << "</g></svg>\n";
}
