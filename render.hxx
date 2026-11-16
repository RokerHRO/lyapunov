#ifndef RENDER_HXX
#define RENDER_HXX

#include "render.hh"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <omp.h>

#include "rgb.hh"

#define CLAMP(x) static_cast<uint8_t>(((x) > 255) ? (255u) : ((x) < 0) ? (0u) : (x))

template<class PXL>
PXL lambda2pxl(double lambda, const FractalParams& params);

template<>
RGB lambda2pxl<RGB>(const double lambda, const FractalParams& params)
{
	double r=0,g=0,b=0;
	if (lambda > 0)
	{
		b = lambda/params.lambda_max;
	} else {
		r = 1 - pow(lambda/params.lambda_min, 2/3.0);
		g = 1 - pow(lambda/params.lambda_min, 1/3.0);
		b = 0;
	}

	return RGB{
		CLAMP(r*255),
		CLAMP(g*255),
		CLAMP(b*255)
	};
}

template<>
uint8_t lambda2pxl<uint8_t>(const double lambda, const FractalParams& /*params*/)
{
	return CLAMP( 127 - lambda*100 );
}


// extension: 3rd dimension: value c is fixed, selected by 'C' in sequence string
template<class PXL>
void render_lyapunov(std::vector<PXL>& img, const FractalParams& params, const FractalType& ftype, const unsigned frame)
{
	const int height = int(params.a.size);
	const int width  = int(params.b.size);

	constexpr unsigned THREADS = 32;
	double tvals[THREADS][5] = { };

	for (auto & tval : tvals)
	{
		ftype.per_image(tval, frame);
	}

#pragma omp parallel for schedule(dynamic)
	for (int ai = 0; ai <height; ++ai)
	{
		const int tn = omp_get_thread_num();
		double* vals = tvals[tn];
		ftype.per_line(vals, frame, ai);

		fprintf(stderr, "\r%4d of %4d   ", ai, height);
		for (int bi = 0; bi < width; ++bi)
		{
			ftype.per_pixel(vals, frame, ai, bi);

			double x = 0.5;
			/* one round without derivating, so that the value 0.5 is avoided */
			for (unsigned m = 0; m < params.sequence.length(); m++)
			{
				const char ch = params.sequence[m];
				const double r = vals[ch-'A']; //ch=='B' ? b : (ch=='A' ? a : c);
				x = r*x*(1-x);
			}
			double sum_log_deriv = 0;
			double prod_deriv = 1;
			bool diverged = false;
			for (int n = 0; n < params.max_iter; n++)
			{
				/* [0,1] is invariant for r in (0,4): once the orbit leaves it, it
				   runs off to -infinity and lambda -> +infinity (pure blue), so
				   the remaining iterations can be skipped */
				if (x < 0.0 || x > 1.0)
				{
					diverged = true;
					break;
				}
				double round_deriv = 1;
				for (unsigned m = 0; m < params.sequence.length(); m++)
				{
					const char ch = params.sequence[m];
					const double r = vals[ch-'A']; // ch=='B' ? b : (ch=='A' ? a : c);
					round_deriv *= r*(1-2*x);
					x = r*x*(1-x);
				}
				/* accumulate the derivative product over several rounds and take
				   one logarithm only when it drifts out of a comfortable range
				   (one per round is what the old code did, one per pixel usually) */
				prod_deriv *= round_deriv;
				const double a = fabs(prod_deriv);
				if (a > 1e100 || a < 1e-100 || n == params.max_iter-1)
				{
					sum_log_deriv += log(a);
					prod_deriv = 1;
				}
			}
			const double lambda = diverged ? 1e30 : sum_log_deriv / (params.max_iter*params.sequence.length());
			img[size_t(ai) * width + bi] = lambda2pxl<PXL>(lambda, params);
		}
	}

	fputs(" \n", stderr);
}

#endif // RENDER_HXX
