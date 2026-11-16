
#ifndef FRACTAL_PARAMS_HH
#define FRACTAL_PARAMS_HH

#include <string>

// one dimension of the fractal computation: the number of pixel and the
// min/max parameter value of that dimension
struct Dimension
{
	unsigned size = 0; // in pixel
	double min = 0.0;
	double max = 0.0;
};

// all parameters of one fractal computation (shared by the renderer, the
// grid overlay and the file writers); must outlive the FractalType created
// from it (FractalType holds references to the dimensions)
struct FractalParams
{
	Dimension a{ 800, 3.999, 2.4 }; // runs along the lines (y axis)
	Dimension b{ 800, 2.4, 3.999 }; // runs along the columns (x axis)
	Dimension c{ 0, 2.5, 4.0 };     // animation frames / 3rd dimension (3df output); size 0: none

	int max_iter = 1000; // number of rounds
	std::string sequence = "BBBABBAAAAAA"; // each character selects the dimension (A..E) of that logistic map round
	double d = 3.7; // 4th dimension (fixed)
	double e = 3.0; // 5th dimension (fixed)

	// lambda range of the RGB color map; somewhat empirical, in order to
	// match Wickerprint's original colors
	double lambda_min = -2.55;
	double lambda_max = 0.3959;
};

#endif
