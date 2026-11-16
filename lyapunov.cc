#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <cctype>
#include <cmath>
#include <cstring>
#include <string>
#include <getopt.h>
#include <utility>
#include <vector>
#include <omp.h>

#include "fractal_type.hh"
#include "string_case.hh"
#include "rgb.hh"
#include "draw.hh"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#define CLAMP(x) static_cast<uint8_t>(((x) > 255) ? (255u) : ((x) < 0) ? (0u) : (x))

//int seq_length = 12;
//int seq[12] = {1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0};

std::string center;
std::string type = "plain";
std::string format = "pnm";
std::string outfile;
std::string post_frame_cmd;
const char* progname = "lyapunov";

double valueF = 3.0;

bool setC = false;
bool setD = false;
bool setE = false;
bool setF = false;

// all parameters of the fractal computation (documented and defaulted in
// fractal_params.hh)
FractalParams fractal_params;
int gstart = 0; // first frame to render (0-based, of the total frame count)
int gend = 0;   // one past the last frame to render (0: render all frames)

std::vector<double> params;
std::string warp_name;
std::vector<double> warp_params;
bool warp_animated = false;

bool grid = false;
std::string grid_spec;
enum GridSize { GridDense, GridMedium, GridSparse };
const Font* grid_font = nullptr; // set in main()
GridSize grid_size = GridMedium;

FILE* of = nullptr;
std::unique_ptr<FractalType> fractal_type;


RGB grid_fg{ 0xFF, 0xFF, 0xFF };
RGB grid_bg{ 0x00, 0x00, 0x00 };


void run_post_frame_cmd(const std::string& out_name, const unsigned frame_nr)
{
	// {o} is substituted last: the file name itself is already fully expanded
	const std::string cmd = replace_all(fractal_params.expand_patterns(post_frame_cmd, frame_nr), "{o}", out_name);
	fprintf(stderr, "=== Post frame command: %s ===\n", cmd.c_str());
	const int rc = std::system(cmd.c_str());
	if (rc != 0)
	{
		fprintf(stderr, "Post frame command for frame #%u failed (status %d), continuing.\n", frame_nr, rc);
	}
}

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
void lyapunov(std::vector<PXL>& img, const FractalParams& params, const FractalType& ftype, const unsigned frame)
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


// round the raw step up to a "nice" step (1, 2 or 5 times a power of ten)
double grid_step(const double raw)
{
	const double p = pow(10.0, floor(log10(raw)));
	for (const double m : {1.0, 2.0, 5.0})
		if (raw <= m * p)
			return m * p;
	return 10.0 * p;
}

std::string grid_label(const double value, const int decimals)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%.*f", decimals, value);
	return buf;
}

// a or b as a linear function of the continuous pixel position (px, py):
// value = c0 + gx*px + gy*py
struct GridField
{
	double c0, gx, gy;

	double value_at(const double px, const double py) const { return c0 + gx*px + gy*py; }

	// the field's minimum and maximum over the window: a linear function
	// over a rectangle achieves its extrema at a corner, so the four
	// corners are evaluated
	std::pair<double, double> range(const double pxlo, const double pxhi, const double pylo, const double pyhi) const
	{
		double lo = value_at(pxlo, pylo), hi = lo;
		for (const auto& p : {std::pair{pxhi, pylo}, std::pair{pxlo, pyhi}, std::pair{pxhi, pyhi}})
		{
			const double v = value_at(p.first, p.second);
			lo = std::min(lo, v);
			hi = std::max(hi, v);
		}
		return {lo, hi};
	}
};

// intersect the line field(px, py) = value with the window rectangle
bool grid_line_endpoints(const GridField& field, const double value,
		const double pxlo, const double pxhi, const double pylo, const double pyhi,
		double& x0, double& y0, double& x1, double& y1)
{
	const double v = value - field.c0;
	constexpr double EPS = 1e-12;
	double pts[2][2];
	int n = 0;
	if (n < 2 && std::abs(field.gy) > EPS)
		for (const double px : {pxlo, pxhi})
		{
			const double py = (v - field.gx*px) / field.gy;
			if (py >= pylo && py <= pyhi)
			{
				pts[n][0] = px;
				pts[n][1] = py;
				++n;
			}
		}
	if (n < 2 && std::abs(field.gx) > EPS)
		for (const double py : {pylo, pyhi})
		{
			const double px = (v - field.gy*py) / field.gx;
			if (px >= pxlo && px <= pxhi)
			{
				pts[n][0] = px;
				pts[n][1] = py;
				++n;
			}
		}
	if (n < 2)
		return false;
	x0 = pts[0][0]; y0 = pts[0][1];
	x1 = pts[1][0]; y1 = pts[1][1];
	return true;
}

// where to place the label: at the end of the line at the left/bottom/top edge
enum GridLabelPosition { LabelLeft, LabelBottom, LabelTop };

// the parameter axis the tileable grid lines belong to (a: sine over the
// lines, horizontal; b: cosine over the columns, vertical)
enum GridAxis { AxisA, AxisB };

// draw one grid line; if text is non-null, also draw its label at the
// left/bottom/top end of the line
void draw_grid_line_labeled(std::vector<RGB>& img, const int W, const int H,
		const double x0, const double y0, const double x1, const double y1,
		const std::string* text, const GridLabelPosition label_position, const RGB& fg)
{
	draw_line(img, W, H, int(x0), int(y0), int(x1), int(y1), fg);
	if (!text)
		return;
	const int label_width = int(text->size()) * int(grid_font->width);
	const int label_height = int(grid_font->height);
	int lx, ly;
	if (label_position == LabelLeft) // constant-a line: label near the left end
	{
		const bool left = x0 <= x1;
		lx = int(left ? x0 : x1) + 2;
		ly = int(left ? y0 : y1) - label_height / 2;
	}
	else if (label_position == LabelBottom) // constant-b line: label near the bottom end
	{
		const bool bottom = y0 >= y1;
		lx = int(bottom ? x0 : x1) - label_width / 2;
		ly = int(bottom ? y0 : y1) - label_height - 2;
	}
	else // constant-c line: label near the top end
	{
		const bool top = y0 <= y1;
		lx = int(top ? x0 : x1) + 2;
		ly = int(top ? y0 : y1) + 1;
	}
	lx = std::clamp(lx, 1, W - label_width - 1);
	ly = std::clamp(ly, 1, H - label_height - 1);
	draw_shaded_text(img, W, H, lx, ly, *text, fg, grid_bg, *grid_font);
}

// draw coordinate lines at "nice" values of the given parameter (a, b or c),
// with labels; the lines are straight for all supported types
void draw_grid_value_lines(std::vector<RGB>& img, const int W, const int H, const GridField& f,
		const double lo, const double hi, const double step, const int decimals,
		const GridLabelPosition label_position, const RGB& fg)
{
	const double pxlo = -0.5, pxhi = W - 0.5, pylo = -0.5, pyhi = H - 0.5;
	const int n_first = int(ceil(lo / step - 1e-9));
	const int n_last = int(floor(hi / step + 1e-9));

	// label every label_stride-th line: the smallest number for which
	// neighboring labels do not overlap (measured along the edge the labels run on)
	int label_stride = 1;
	const int n_lines = n_last - n_first + 1;
	if (label_position == LabelLeft) // a labels: stacked vertically, must fit their height
	{
		const int label_height = int(grid_font->height);
		label_stride = std::max(1, (label_height * n_lines + H - 1) / H);
	}
	else // b/c labels: run horizontally along the bottom/top edge
	{
		const std::string first_label = grid_label(n_first * step, decimals);
		const std::string last_label = grid_label(n_last * step, decimals);
		const int label_width = int(std::max(first_label.size(), last_label.size())) * int(grid_font->width);
		label_stride = std::max(1, (label_width * n_lines + W - 1) / W);
	}

	for (int n = n_first; n <= n_last; ++n)
	{
		const double value = n * step;
		double x0, y0, x1, y1;
		if (!grid_line_endpoints(f, value, pxlo, pxhi, pylo, pyhi, x0, y0, x1, y1))
			continue;
		const std::string text = grid_label(value, decimals);
		draw_grid_line_labeled(img, W, H, x0, y0, x1, y1,
			(n - n_first) % label_stride ? nullptr : &text, label_position, fg);
	}
}

// tileable / tileable-loop: the window is one full period of a sine (a along
// the lines) and a cosine (b along the columns), so a given value occurs on
// up to two lines; the positions are solved from sin/cos = target
void draw_grid_sinus_axis(std::vector<RGB>& img, const int W, const int H, const GridAxis axis,
		const double center, const double radius, const double lo, const double hi,
		const double step, const int decimals, const GridLabelPosition label_position, const RGB& fg)
{
	// the line positions run in pixel units along the varying axis (py for
	// the a lines, px for the b lines); the pixel window is [-0.5, N-0.5]
	const double pos_lo = -0.5, pos_hi = (axis == AxisA ? H : W) - 0.5;
	// first/last line index: the grid values are n * step for n_first..n_last
	const int n_first = int(ceil(lo / step - 1e-9));
	const int n_last = int(floor(hi / step + 1e-9));

	// label every label_stride-th value so that neighboring labels do not
	// overlap along the labeled edge; a value can have two lines, so double
	// the line count
	const int label_stride = [&] {
		const int n_lines = 2 * (n_last - n_first + 1);
		if (label_position == LabelLeft)
		{
			const int label_height = int(grid_font->height);
			return std::max(1, (label_height * n_lines + H - 1) / H);
		}
		const std::string first_label = grid_label(n_first * step, decimals);
		const std::string last_label = grid_label(n_last * step, decimals);
		const int label_width = int(std::max(first_label.size(), last_label.size())) * int(grid_font->width);
		return std::max(1, (label_width * n_lines + W - 1) / W);
	}();

	for (int n = n_first; n <= n_last; ++n)
	{
		const double value = n * step;
		// the sine/cosine value the window function must take for this line
		const double target = (value - center) / radius;
		if (target < -1.0 - 1e-9 || target > 1.0 + 1e-9)
			continue;
		const double target_clamped = std::clamp(target, -1.0, 1.0);
		// the two angle solutions of sin(angle) = target / cos(angle) = target
		// in one period; the first is wrapped by 2*pi if it falls below the
		// window start
		const double scale = axis == AxisA ? H : W;
		double angles[2] = { axis == AxisA ? asin(target_clamped) : acos(target_clamped), axis == AxisA ? M_PI - asin(target_clamped) : 2*M_PI - acos(target_clamped) };
		if (angles[0] * scale / (2*M_PI) < pos_lo)
			angles[0] += 2*M_PI;
		if (angles[0] > angles[1])
			std::swap(angles[0], angles[1]);
		const std::string text = grid_label(value, decimals);
		double last = -1e30;
		for (const double angle : angles)
		{
			const double pos = angle * scale / (2*M_PI);
			if (pos < pos_lo || pos > pos_hi || pos - last < 1e-9) // the solutions coincide at s = +/-1
				continue;
			last = pos;
			double x0, y0, x1, y1;
			if (axis == AxisA) // horizontal line
			{
				x0 = -0.5; y0 = pos; x1 = W - 0.5; y1 = pos;
			}
			else // vertical line
			{
				x0 = pos; y0 = -0.5; x1 = pos; y1 = H - 0.5;
			}
			draw_grid_line_labeled(img, W, H, x0, y0, x1, y1,
				(n - n_first) % label_stride == 0 ? &text : nullptr, label_position, fg);
		}
	}
}

void draw_grid(std::vector<RGB>& img)
{
	// the grid is only supported for types with a straight-line mapping from
	// the pixel position to the (a, b, c) parameter space: affine for
	// plain/circle/zoom/flights, one sine/cosine period per axis for tileable
	switch (case_hash(fractal_type->name()))
	{
		case "plain"_case:  [[fallthrough]];
		case "circle"_case: [[fallthrough]];
		case "zoom"_case:   [[fallthrough]];
		case "tileable"_case: [[fallthrough]];
		case "tileable-loop"_case: [[fallthrough]];
		case "six_side"_case: [[fallthrough]];
		case "wander"_case: [[fallthrough]];
		case "rotor"_case:  break;
		default:
			return;
	}
	// a warp distorts u, v before the (affine) mapping: the lines would be curved
	if (!warp_name.empty())
		return;

	const int W = int(fractal_params.b.size), H = int(fractal_params.a.size);
	const double pxlo = -0.5, pxhi = W - 0.5, pylo = -0.5, pyhi = H - 0.5;

	bool flight = false, tileable = false;
	switch (case_hash(fractal_type->name()))
	{
		case "six_side"_case: [[fallthrough]];
		case "wander"_case:   [[fallthrough]];
		case "rotor"_case:    flight = true; break;
		case "tileable"_case: [[fallthrough]];
		case "tileable-loop"_case: tileable = true; break;
		default: break;
	}

	// (a, b, c) at the continuous pixel position (px, py); c is only defined
	// for the flight types (for the other types it is constant per frame)
	const auto abc_at = [&](const double px, const double py, double& a, double& b, double& c)
	{
		switch (case_hash(fractal_type->name()))
		{
			case "zoom"_case:
			{
				// the same (x, y) as in Zoom::per_pixel, continuous in px, py
				const auto& z = static_cast<const Zoom&>(*fractal_type);
				z.window_ab((px + 0.5) / W - 0.5, (py + 0.5) / H - 0.5, a, b);
				break;
			}
			case "six_side"_case: [[fallthrough]];
			case "wander"_case: [[fallthrough]];
			case "rotor"_case:
			{
				// the same (u, v) as in Flight::per_pixel (without warp)
				const unsigned side = H > W ? H : W;
				const double u = (px + 0.5 + (side - W) / 2) / side;
				const double v = (py + 0.5 + (side - H) / 2) / side;
				const auto& f = static_cast<const Flight&>(*fractal_type);
				f.window_ab(u, v, a, b);
				c = f.window_c(u, v);
				break;
			}
			default:
				// plain / circle: a runs along the lines, b along the columns
				// (the same (line, column) mapping as in SimpleFractalType)
				a = fractal_params.a.min + (fractal_params.a.max-fractal_params.a.min) * (py + 0.5) / H;
				b = fractal_params.b.min + (fractal_params.b.max-fractal_params.b.min) * (px + 0.5) / W;
				break;
		}
		if (!flight)
			c = 0.0;
	};

	// linear coefficients by sampling: c0 at (0, 0), gradient from (1, 0) and (0, 1)
	auto make_field = [&](const int which) -> GridField
	{
		double v0[3], v1[3], v2[3];
		abc_at(0, 0, v0[0], v0[1], v0[2]);
		abc_at(1, 0, v1[0], v1[1], v1[2]);
		abc_at(0, 1, v2[0], v2[1], v2[2]);
		return {v0[which], v1[which] - v0[which], v2[which] - v0[which]};
	};
	const GridField fa = tileable ? GridField{} : make_field(0);
	const GridField fb = tileable ? GridField{} : make_field(1);
	const GridField fc = flight ? make_field(2) : GridField{};

	std::pair<double, double> ra, rb;
	if (tileable) // one full period: the window covers the whole a/b range
	{
		ra = { std::min(fractal_params.a.min, fractal_params.a.max), std::max(fractal_params.a.min, fractal_params.a.max) };
		rb = { std::min(fractal_params.b.min, fractal_params.b.max), std::max(fractal_params.b.min, fractal_params.b.max) };
	}
	else
	{
		ra = fa.range(pxlo, pxhi, pylo, pyhi);
		rb = fb.range(pxlo, pxhi, pylo, pyhi);
	}
	const auto rc = fc.range(pxlo, pxhi, pylo, pyhi);

	// target number of lines on the larger axis (grid_size)
	const int target = grid_size == GridDense ? 40 : grid_size == GridSparse ? 8 : 12;
	const double span = std::max(ra.second - ra.first, rb.second - rb.first);
	const double step = grid_step(span / target);
	const int decimals = std::max(0, int(ceil(-log10(step) - 1e-9)));

	// per-axis colors derived from grid_fg: b halves the red channel, c halves
	// the green channel (white fg -> #7fffff and #ff7fff)
	const RGB fg_b{ uint8_t(grid_fg.r/2), grid_fg.g, grid_fg.b };
	const RGB fg_c{ grid_fg.r, uint8_t(grid_fg.g/2), grid_fg.b };

	if (flight)
		draw_grid_value_lines(img, W, H, fc, rc.first, rc.second, step, decimals, LabelTop, fg_c);
	if (tileable)
	{
		draw_grid_sinus_axis(img, W, H, AxisB, (fractal_params.b.min + fractal_params.b.max) / 2, (fractal_params.b.max - fractal_params.b.min) / 2, rb.first, rb.second, step, decimals, LabelBottom, fg_b);
		draw_grid_sinus_axis(img, W, H, AxisA, (fractal_params.a.min + fractal_params.a.max) / 2, (fractal_params.a.max - fractal_params.a.min) / 2, ra.first, ra.second, step, decimals, LabelLeft, grid_fg);
	}
	else
	{
		draw_grid_value_lines(img, W, H, fb, rb.first, rb.second, step, decimals, LabelBottom, fg_b);
		draw_grid_value_lines(img, W, H, fa, ra.first, ra.second, step, decimals, LabelLeft, grid_fg);
	}
}

void make_ppm(const unsigned frame_nr)
{
	std::string outfile_ppm;
	if (!outfile.empty())
	{
		outfile_ppm = fractal_params.expand_patterns(outfile, frame_nr);
		if (outfile_ppm.empty() || outfile_ppm.back() == '/' || outfile_ppm.back() == '\\')
			throw std::runtime_error("the output filename is empty after replacing the template patterns");
		if (!outfile_ppm.ends_with(".ppm"))
			outfile_ppm += ".ppm";

		of = fopen(outfile_ppm.c_str(), "wb");
		if (!of)
		{
			throw std::runtime_error("Cannot open output file \"" + outfile_ppm + "\"");
		}

		fprintf(stderr, "=== Write frame #%u to file \"%s\" ===\n", frame_nr, outfile_ppm.c_str());
	}else
	{
		fprintf(stderr, "=== Write frame #%u to stdout ===\n", frame_nr);
	}

	// RGB image:
	std::vector<RGB> img;
	img.resize(fractal_params.a.size * fractal_params.b.size);

	// PPM header
	fprintf( of,"P6\n");
	fflush(of);
	lyapunov(img, fractal_params, *fractal_type, frame_nr);
	if (grid)
		draw_grid(img);
	fprintf(of, "# Lyapunov: '%s' max_iter=%d  seq='%s'  a=%f ... %f  b=%f ... %f ",
		fractal_type->name(),
		fractal_params.max_iter, fractal_params.sequence.c_str(),
		fractal_params.a.min, fractal_params.a.max, fractal_params.b.min, fractal_params.b.max);

	if(setC) fprintf(of, " frame=%u/%u  C=%f ... %f", frame_nr, fractal_params.c.size, fractal_params.c.min, fractal_params.c.max);
	if(setD) fprintf(of, " D=%f", fractal_params.d);
	if(setE) fprintf(of, " E=%f", fractal_params.e);
	if(setF) fprintf(of, " F=%f", valueF);
	
	fprintf(of, "\n"
		"%d %d 255\n",
		int(fractal_params.b.size), int(fractal_params.a.size)
		);
	
	for (int ai = 0; ai < int(fractal_params.a.size); ++ai)
	{
		fwrite(&img[size_t(ai) * fractal_params.b.size], sizeof(RGB), fractal_params.b.size, of);
	}
	fflush(of);

	if (!outfile.empty())
	{
		fclose(of);
		if (!post_frame_cmd.empty())
			run_post_frame_cmd(outfile_ppm, frame_nr);
	}
}

void make_ppm()
{
	if (fractal_params.c.size==0)
	{
		make_ppm(0);
	}else for (int frame_nr=gstart; frame_nr < gend; ++frame_nr)
	{
		make_ppm(frame_nr);
	}
}

inline
uint8_t get(const std::vector<std::vector<uint8_t>>& img, const int x, const int y, const int z)
{
	if( unsigned(y) >= fractal_params.a.size ) return 0;
	if( unsigned(x) >= fractal_params.b.size ) return 0;
	return img[z][size_t(y) * fractal_params.b.size + x];
}


void make_3df()
{
	if (!outfile.empty())
	{
		std::string outfile_3df = outfile;
		if (!outfile_3df.ends_with(".3df"))
			outfile_3df += ".3df";
		of = fopen(outfile_3df.c_str(), "wb");
		if (!of)
		{
			throw std::runtime_error("Cannot open output file \"" + outfile_3df + "\"");
		}
	}
	// 3df header
	fprintf(of, "%c%c%c%c%c%c",
		fractal_params.b.size >> 8, fractal_params.b.size & 255,
		fractal_params.a.size >> 8, fractal_params.a.size & 255,
		fractal_params.c.size >> 8, fractal_params.c.size & 255
	);
	fflush(of);
	
	// 3 images:
	std::vector<std::vector<uint8_t>> img;
	img.resize(3, std::vector<uint8_t>(fractal_params.a.size * fractal_params.b.size));
	
	for(int ci=0; ci<int(fractal_params.c.size); ++ci)
	{
		const double c = (fractal_params.c.size > 1) ? fractal_params.c.min + (fractal_params.c.max-fractal_params.c.min)/int(fractal_params.c.size-1)*ci : (fractal_params.c.min+fractal_params.c.max)/2;
		fprintf(stderr, "Frame %i of %i. c=%f \n", ci, int(fractal_params.c.size), c);
		lyapunov(img[2], fractal_params, *fractal_type, ci);
		
		fprintf(stderr,"\tblurring...\n");
		for(int y=0; y<int(fractal_params.a.size); ++y)
		for(int x=0; x<int(fractal_params.b.size); ++x)
		{
			unsigned u = 0;
			for(int z=0; z<3; ++z)
			for(int dy=-1; dy<=+1; ++dy)
			for(int dx=-1; dx<=+1; ++dx)
			{
				u += get(img, x+dx, y+dy, z);
			}
			const uint8_t u8 = u/27;
			fputc(u8, of);
		}
		fflush(of);
		img[0].swap(img[1]);
		img[1].swap(img[2]);
	}

	if (!outfile.empty())
	{
		fclose(of);
	}
}

void print_type_help(FILE* f)
{
	fprintf(f, "Fractal types:\n"
		"\tplain  : A=y axis, B=x axis, C=z axis (animation or 3rd dimension), D=fix, E=fix\n"
		"\t6sides : tilts the (A,B) slice at C=cmin about 6 cube edges (90 degrees each,\n"
		"\t         axes A,C,B,A,C,B), visiting all 6 faces of the ABC cube. D=fix, E=fix\n"
		"\t         -p [rotate][,ease]: rotate 1: the window rotates 180 degrees in the\n"
		"\t             plane over the animation (seamless loop); ease 1: smoothstep\n"
		"\t             timing, the window settles at every face; ease 2: quaternion\n"
		"\t             Catmull-Rom path through the faces (continuous axis blend,\n"
		"\t             cruise, no pause at the faces) (default 0,0)\n"
		"\twander : the window center wanders along a closed Catmull-Rom spline through N random\n"
		"\t         keyframes (random orientations, constant cruise speed), seamless loop\n"
		"\t         -p [keyframes][,seed][,scale]  (default 6,0,1.0; keyframes >= 3, scale = window\n"
		"\t             size relative to a cube face, < 1 = zoom)\n"
		"\trotor  : the (A,B) slice at C=0.5 (A vertical, B horizontal) spins in place about an axis\n"
		"\t         through the cube center, whole number of full turns over the animation, seamless loop\n"
		"\t         -p [ax,ay,az][,turns]  (default 1,0,0,1: the A axis, one turn; e.g. 1,1,1 = the\n"
		"\t             space diagonal)\n"
		"\tcircle : C+D circle\n"
		"\tsphere : Sphere in the ABC cube. D=fix, E=fix\n"
		"\ttileable : A and B each trace a full circle, so the image tiles seamlessly in both directions\n"
		"\ttileable-loop : tileable, plus C+D circle over the frames, so the animation loops seamlessly\n"
		"\tzoom : fly-through into the (A,B) plane: over the whole animation the rectangle is scaled about its center and rotated\n"
		"\t            -p <total_zoom_factor>,<total_rotation_degrees>  (default 0.1,360; 0<f<1: zoom in, f>1: zoom out)\n"
		"\twavy : plain (A,B) plane with a sinusoidal warp: A wobbles along the columns, B along the lines\n"
		"\t            -p <freq_a>,<amp_a>,<phase_a_deg>[,<freq_b>,<amp_b>,<phase_b_deg>]  (frequency = waves per axis,\n"
		"\t                 phase in degrees, 360 = one wave; B triple defaults to the A values; default 1,0.02,0)\n"
		"\tnoise : plain (A,B) plane warped by 2D Perlin noise\n"
		"\t            -p <frequency>,<amp_a>,<amp_b>[,<seed>]  (frequency = noise cells per axis, amplitude =\n"
		"\t                 fraction of image size; default 3,0.333,0.333,0)\n"
		"\tinversion : plain (A,B) plane, but the points are reflected about a fixed circle (P' = C + R^2*(P-C)/|P-C|^2)\n"
		"\t            -p [center_a,center_b[,radius]]  (default: window center, half of the A range)\n"
		"\n"
		"Warps (combine with a flight type via -w, parameters via -P):\n"
		"\ta distortion of the image window, applied before the flight mapping, so it\n"
		"\tcan be combined with 6sides, wander and rotor:\n"
		"\t-t <flight> -w <warp> -p <flight params> -P <warp params> [-A]\n"
		"\t-A        : animate the warp: its main effect (amplitude, twist angle, lens\n"
		"\t             center) ramps up linearly from 0 (first frame) to the full -P value\n"
		"\t             (last frame); supported by wavy, noise, mobius, swirl and twist\n"
		"\t             (zoom is already animated, inversion has no sensible ramp).\n"
		"\t             For wavy the phases instead sweep one full wave (360 degrees):\n"
		"\t             the wave travels, seamless loop for integer frequencies\n"
		"\twavy      : sinusoidal warp (same parameters as the wavy type, amplitudes as a fraction of the window)\n"
		"\tnoise     : 2D Perlin noise warp (same parameters as the noise type, amplitudes as a fraction of the window)\n"
		"\tinversion : reflection of the window about a fixed circle (center and radius as a fraction of the window,\n"
		"\t             default 0.5,0.5,0.5)\n"
		"\tmobius    : hyperbolic lens (Poincare-disk Moebius transform): pulls the content toward (tx,ty),\n"
		"\t             the opposite side is compressed conformally -P [tx,ty][,radius]  (default 0.15,0,0.5)\n"
		"\tswirl     : vortex about the window center, twist angle grows exponentially with the distance\n"
		"\t             (center almost undisturbed, border strongly wound up)\n"
		"\t             -P <twist_degrees>[,radius][,decay]  (default 180,0.7071,3)\n"
		"\ttwist     : rotate the window about its center, angle proportional to the square of the distance\n"
		"\t             (center fixed, border winds up) -P <twist_degrees>[,radius]  (default 180,0.7071)\n"
		"\tzoom      : scale and rotate the window about its center over the animation (same parameters as the zoom type)\n"
	);
}

void print_format_help(FILE* f)
{
	fprintf(f, "Output formats:\n"
		"\tppm : a PPM file or series of PPM files (if -o is given) or concatenated series of PPM images to stdout (of no -o given)\n"
		"\t3df : a 3DF voxel file\n"
	);
}

std::string font_names()
{
	std::string list;
	for (const auto& f : fonts)
	{
		if (!list.empty()) list += ", ";
		list += f.first;
	}
	return list;
}

void print_help(FILE* f)
{
	std::string font_list = font_names();
	font_list += ", default: " + fonts.cbegin()->first;

	fprintf(f,
		"Generate a Lyapunov fractal as PPM or 3df to stdout.\n\n"
		"Usage: %s [-W width] [-H height] [-f frames] [-s sequence] [-i max_iter]\n"
		"     [-t fractal_type ] [-o output_filename_base] [-O format]\n"
		"   { [-c center_x:center_y:size] | [-x min_x] [-X max_x] [-y min_y] [-Y max_y] } [-z min_z] [-Z max_z]\n"
		"     [-C value] [-D value] [-E value] [-F value]\n"
		"\n"
		"\t-h, --help : print this help and exit\n"
		"\t-W width, --width : width of image in pixel (default: %i)\n"
		"\t-H height, --height : height of image in pixel (default: %i)\n"
		"\t-f frames, --frames : number of frames or size in in z direction\n"
		"\t-G start:end, --frame-range : render only the frames in the range [start, end) of the -f frames\n"
		"\t                               (frame numbers and timing use the total -f count, useful for\n"
		"\t                               rendering chunks in parallel; default: all frames)\n"
		"\t-s seq, --sequence : AB or ABC sequence for fractal (default: \"%s\")\n"
		"\t-i iter, --iterations : iterations per pixel (default: %i)\n"
		"\n"
		"\t-t type, --type : fractal type (default: %s, run \"-t help\" to get a list of types)\n"
		"\t-o outfile, --output : basename of the output file (stdout if not given)\n"
		"\t                      These patterns are replaced: {s}: sequence, {i}: iter, {f}: frame_nr\n"
		"\t-r command, --post-frame-cmd : run a shell command after each written frame (pnm/ppm only, requires -o)\n"
		"\t                               Patterns: {s}: sequence, {i}: iter, {f}: frame_nr, {o}: output file name\n"
		"\t-O format, --format : either \"pnm\" for image (series) or \"3df\" for 3D voxel file (default: \"%s\")\n"
		"\n"
		"\t-c cx:cy:x_size, --center : center coordinates and size.\n"
		"\t   The x, X, y and Y values will be calculated from them automatically:\n"
		"\t-x min_x, --min-x : minimum x value (default: %f)\n"
		"\t-X max_x, --max-x : maximum x value (default: %f)\n"
		"\t-y min_y, --min-y : minimum y value (default: %f)\n"
		"\t-Y max_y, --max-y : maximum y value (default: %f)\n"
		"\t-z min_z, --min-z : 3rd dimension: (start) value of z. (selected by \"C\" in sequence string. default: %f)\n"
		"\t-Z max_z, --max-z : maximum value of z (default: %f) \n"
		"\t-C value, --dim-c : value for 3th dimension (default %f)\n"
		"\t-D value, --dim-d : value for 4th dimension (default %f)\n"
		"\t-E value, --dim-e : value for 5th dimension (default %f)\n"
		"\n"
		"\t-p params[,...], --params : comma-separated list of additional numerical parameters for -t\n"
		"\t-w warp, --warp : warp of the image window (only with the flight types 6sides, wander, rotor;\n"
		"\t                  supported: wavy, noise, inversion, mobius, swirl, twist, zoom)\n"
		"\t-P params[,...], --warp-params : comma-separated list of numerical parameters for -w\n"
		"\t-A, --animate-warp : animate the warp over the frames (see -t help)\n"
		"\t--grid[=font,size,fg,bg] : draw and label coordinate lines on the rendered PPM frames\n"
		"\t                           (supported: plain, circle, zoom, tileable, tileable-loop, and the flight types\n"
		"\t                           6sides, wander, rotor; font: %s; grid size dense, medium or sparse, default medium;\n"
		"\t                           fg/bg #RRGGBB, default #FFFFFF/#000000; trailing parameters may be omitted).\n"
		"\t                           Colors: a = fg, b = fg with the red channel halved (#7fffff for white);\n"
		"\t                           for the flight types also c lines are drawn, using fg with the green\n"
		"\t                           channel halved (#ff7fff for white)\n"
		, progname, int(fractal_params.b.size), int(fractal_params.a.size), fractal_params.sequence.c_str(), fractal_params.max_iter,
		type.c_str(), format.c_str(),
		fractal_params.b.min, fractal_params.b.max, fractal_params.a.min, fractal_params.a.max,
		fractal_params.c.min, fractal_params.c.max, fractal_params.c.min, fractal_params.d, fractal_params.e,
		font_list.c_str()
		);
}

std::vector<double> parse_params(const char* p)
{
	std::vector<double> params;
	std::istringstream ss(p);
	double value;
	while (ss >> value)
	{
		params.push_back(value);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return params;
}

RGB parse_hex_color(const std::string& t)
{
	if (t.size() != 7 || t[0] != '#')
		throw std::invalid_argument("Invalid grid color \"" + t + "\" (expected #RRGGBB).");
	for (size_t i = 1; i < t.size(); ++i)
		if (!isxdigit(uint8_t(t[i])))
			throw std::invalid_argument("Invalid grid color \"" + t + "\" (expected #RRGGBB).");
	return RGB{
		uint8_t(std::stoul(t.substr(1, 2), nullptr, 16)),
		uint8_t(std::stoul(t.substr(3, 2), nullptr, 16)),
		uint8_t(std::stoul(t.substr(5, 2), nullptr, 16))
	};
}

// grid spec: [font,size,fg,bg], all parameters optional (in this order), e.g. "24,sparse"
void parse_grid_spec(const std::string& spec)
{
	if (spec.empty())
		return; // getopt_long may pass "" for "--grid" without a value: all defaults
	std::vector<std::string> tokens;
	std::string tok;
	for (char ch : (spec + ','))
	{
		if (ch == ',')
		{
			tokens.push_back(tok);
			tok.clear();
		}
		else
			tok += ch;
	}
	if (tokens.size() > 4)
		throw std::invalid_argument("Grid specification \"" + spec + "\" has too many parameters (max: font,size,fg,bg).");
	for (size_t i = 0; i < tokens.size(); ++i)
	{
		const std::string& t = tokens[i];
		if (t.empty())
			throw std::invalid_argument("Grid parameter " + std::to_string(i+1) + " of \"" + spec + "\" is empty.");
		switch (i)
		{
			case 0:
			{
				const auto f = fonts.find(t);
				if (f == fonts.end())
					throw std::invalid_argument("Invalid grid font \"" + t + "\" (allowed: " + font_names() + ").");
				grid_font = f->second;
				break;
			}
			case 1:
				switch (case_hash(t))
				{
					case "dense"_case:  grid_size = GridDense;  break;
					case "medium"_case: grid_size = GridMedium; break;
					case "sparse"_case: grid_size = GridSparse; break;
					default:
						throw std::invalid_argument("Invalid grid size \"" + t + "\" (allowed: dense, medium, sparse).");
				}
				break;
			case 2: grid_fg = parse_hex_color(t); break;
			case 3: grid_bg = parse_hex_color(t); break;
		}
	}
}


// long-only options (no short form): values >= 256 can never collide with a short option
enum LongOptions { OptGrid = 256 };

static const struct option long_options[] =
{
	{ "width",          required_argument, nullptr, 'W' },
	{ "height",         required_argument, nullptr, 'H' },
	{ "frames",         required_argument, nullptr, 'f' },
	{ "frame-range",    required_argument, nullptr, 'G' },
	{ "sequence",       required_argument, nullptr, 's' },
	{ "iterations",     required_argument, nullptr, 'i' },
	{ "type",           required_argument, nullptr, 't' },
	{ "output",         required_argument, nullptr, 'o' },
	{ "post-frame-cmd", required_argument, nullptr, 'r' },
	{ "format",         required_argument, nullptr, 'O' },
	{ "center",         required_argument, nullptr, 'c' },
	{ "min-x",          required_argument, nullptr, 'x' },
	{ "max-x",          required_argument, nullptr, 'X' },
	{ "min-y",          required_argument, nullptr, 'y' },
	{ "max-y",          required_argument, nullptr, 'Y' },
	{ "min-z",          required_argument, nullptr, 'z' },
	{ "max-z",          required_argument, nullptr, 'Z' },
	{ "dim-c",          required_argument, nullptr, 'C' },
	{ "dim-d",          required_argument, nullptr, 'D' },
	{ "dim-e",          required_argument, nullptr, 'E' },
	{ "dim-f",          required_argument, nullptr, 'F' },
	{ "params",         required_argument, nullptr, 'p' },
	{ "warp",           required_argument, nullptr, 'w' },
	{ "warp-params",    required_argument, nullptr, 'P' },
	{ "animate-warp",   no_argument,       nullptr, 'A' },
	{ "grid",           optional_argument, nullptr, OptGrid },
	{ "help",           no_argument,       nullptr, 'h' },
	{ nullptr,          0,                 nullptr, 0    }
};

int main(int argc, char** argv) try
{
	const auto program_start = std::chrono::system_clock::now();
	progname = argv[0];

	int a;
	while( (a=getopt_long(argc, argv, "hW:H:f:c:x:y:z:X:Y:Z:s:i:o:O:t:C:D:E:F:p:G:r:w:P:A", long_options, nullptr)) != -1 )
	{
		switch(a)
		{
			case 'W':
			{
				const int w = atoi(optarg);
				if (w <= 0)
					throw std::invalid_argument("Option -W/--width expects a positive width in pixel.");
				fractal_params.b.size = unsigned(w);
				break;
			}
			case 'H':
			{
				const int h = atoi(optarg);
				if (h <= 0)
					throw std::invalid_argument("Option -H/--height expects a positive height in pixel.");
				fractal_params.a.size = unsigned(h);
				break;
			}
			case 'f':
			{
				const int f = atoi(optarg);
				if (f < 0)
					throw std::invalid_argument("Option -f/--frames expects a non-negative frame count.");
				fractal_params.c.size = unsigned(f);
				setC = true;
				break;
			}
			case 'z': fractal_params.c.min  = atof(optarg); setC = true; break;
			case 'Z': fractal_params.c.max  = atof(optarg); setC = true; break;
			case 'c': center = optarg;      break;
			case 'x': fractal_params.b.min = atof(optarg);  break;
			case 'y': fractal_params.a.min = atof(optarg);  break;
			case 'X': fractal_params.b.max = atof(optarg);  break;
			case 'Y': fractal_params.a.max = atof(optarg);  break;
			case 's': fractal_params.sequence = optarg;     break;
			case 'i': fractal_params.max_iter = atoi(optarg); break;
			case 't': type = optarg;        break;
			case 'o': outfile = optarg;     break;
			case 'O': format = optarg;      break;
			case 'C': fractal_params.c.min   = atof(optarg); setC = true; break;
			case 'D': fractal_params.d = atof(optarg); setD = true; break;
			case 'E': fractal_params.e = atof(optarg); setE = true; break;
			case 'F': valueF = atof(optarg); setF = true; break;
			case 'p': params = parse_params(optarg); break;
			case 'r': post_frame_cmd = optarg; break;
			case 'w': warp_name = optarg;   break;
			case 'P': warp_params = parse_params(optarg); break;
			case 'A': warp_animated = true;  break;
			case OptGrid: grid = true; if (optarg) grid_spec = optarg; break;
			case 'G':
			{
				const char* colon = strchr(optarg, ':');
				if (colon == nullptr)
					throw std::invalid_argument("Option -G/--frame-range expects start:end.");
				gstart = atoi(optarg);
				gend = atoi(colon+1);
				break;
			}
			case 'h':
				print_help(stdout);
				return 0;
			case '?': [[fallthrough]];
			default:
				print_help(stderr);
				return 1;
		}
	}
	
	if(fractal_params.c.size > 0)
	{
		if(gend == 0) gend = int(fractal_params.c.size); // no -G: render all frames
		if(gstart < 0 || gend <= gstart || gend > int(fractal_params.c.size))
			throw std::invalid_argument("Invalid frame range -G " + std::to_string(gstart) + ":" + std::to_string(gend) +
			                             " (total frames: " + std::to_string(int(fractal_params.c.size)) + ").");
	}

	if(fractal_params.sequence.empty())
		throw std::invalid_argument("Sequence must not be empty.");

	const auto illegal = std::find_if(fractal_params.sequence.begin(), fractal_params.sequence.end(), [](char ch) { return ch < 'A' || ch > 'E'; });
	if (illegal != fractal_params.sequence.end())
		throw std::invalid_argument("Sequence \"" + fractal_params.sequence + "\" contains illegal character '" + *illegal + "' (only 'A'..'E' are allowed).");

	if(fractal_params.max_iter <= 0)
		throw std::invalid_argument("max_iter must be > 0 (got " + std::to_string(fractal_params.max_iter) + ").");

	if(grid)
	{
		grid_font = fonts.cbegin()->second; // default: first entry of fonts
		parse_grid_spec(grid_spec);
	}

	if(!center.empty())
	{
		double cx=0.0, cy=0.0, dx=0.0;
		const int ret = sscanf(center.c_str(), "%lg:%lg:%lg", &cx, &cy, &dx);
		if(ret != 3)
			throw std::invalid_argument("Cannot parse \"" + center + "\" as \"cx:cy:x_size\" parameter.");
		
		const double dy = (dx * int(fractal_params.a.size)) / int(fractal_params.b.size);
		fractal_params.b.min = cx - dx/2;
		fractal_params.b.max = cx + dx/2;
		fractal_params.a.min = cy + dy/2;
		fractal_params.a.max = cy - dy/2;
		
		fprintf(stderr, "Interval calculated from \"-c %s\" ", center.c_str());
	}

	fprintf(stderr, " min_x=%11.9f  max_x=%11.9f  min_y=%11.9f  max_y=%11.9f",
		fractal_params.b.min, fractal_params.b.max, fractal_params.a.min, fractal_params.a.max);

	if(setC) fprintf(stderr, " C=%f ... %f", fractal_params.c.min, fractal_params.c.max);
	if(setD) fprintf(stderr, " D=%f", fractal_params.d);
	if(setE) fprintf(stderr, " E=%f", fractal_params.e);
	if(setF) fprintf(stderr, " F=%f", valueF);
	if (!params.empty())
	{
		fprintf(stderr, " P=%f", params[0]);
		for (unsigned u=1; u<params.size(); ++u)
		{
			fprintf(stderr, ", %f", params[u]);
		}
	}
	if (!warp_name.empty())
	{
		fprintf(stderr, " W=%s", warp_name.c_str());
		if (warp_animated)
		{
			fprintf(stderr, " (animated)");
		}
		if (!warp_params.empty())
		{
			fprintf(stderr, " WP=%f", warp_params[0]);
			for (unsigned u=1; u<warp_params.size(); ++u)
			{
				fprintf(stderr, ", %f", warp_params[u]);
			}
		}
	}
	if (grid)
	{
		fprintf(stderr, " GRID=%s", grid_spec.empty() ? "default" : grid_spec.c_str());
	}

	fputc('\n', stderr);


	fprintf(stderr, "Run with %d threads...\n", omp_get_max_threads());
	of = stdout;

#ifdef _WIN32
	setmode(fileno(stdout),O_BINARY);
	setmode(fileno(stdin),O_BINARY);
#endif

	switch (case_hash(type))
	{
		case "help"_case   : [[fallthrough]];
		case "-h"_case     : [[fallthrough]];
		case "--help"_case : print_type_help(stdout); return 0;
	}

	if (!warp_name.empty())
	{
		switch (case_hash(type))
		{
			case "6sides"_case: [[fallthrough]];
			case "wander"_case: [[fallthrough]];
			case "rotor"_case:  break;
			default:
				throw std::invalid_argument("Warp \"-w " + warp_name + "\" can only be combined with the flight types 6sides, wander and rotor (got type \"" + type + "\").");
		}
	}
	if (warp_animated && warp_name.empty())
		throw std::invalid_argument("Option -A (animate warp) requires -w <warp>.");

	std::unique_ptr<Warp> warp = warp_name.empty() ? nullptr : make_warp(warp_name, warp_params, warp_animated);
	fractal_type = make_fractal_type(type, fractal_params.a, fractal_params.b, fractal_params.c,
		fractal_params.d, fractal_params.e, params, std::move(warp));

	if (case_hash(format) == "3df"_case && fractal_params.c.size == 0)
		throw std::invalid_argument("3df output requires a z dimension: use -f <frames>.");

	if (!post_frame_cmd.empty())
	{
		if (outfile.empty())
			throw std::invalid_argument("Option -r (post frame command) requires -o <outfile>.");
		if (case_hash(format) == "3df"_case)
			throw std::invalid_argument("Option -r (post frame command) is only supported with pnm/ppm output.");
	}

	switch (case_hash(format))
	{
		case "pnm"_case   : [[fallthrough]];
		case "ppm"_case   : make_ppm(); break;

		case "3df"_case   : make_3df(); break;

		case "help"_case  : [[fallthrough]];
		case "-h"_case    : [[fallthrough]];
		case "--help"_case: print_format_help(stdout); return 0;
		default: print_format_help(stderr); return 1;
	}

	fflush(stdout);
	const auto program_end = std::chrono::system_clock::now();
	const std::chrono::hh_mm_ss hms{ program_end - program_start};
	std::stringstream ss;
	ss << hms;

	fprintf(stderr, "All done.\nTotal time: %s \n", ss.str().c_str());
}
catch (const std::invalid_argument& e)
{
	fprintf(stderr, "%s\n", e.what());
	return 2;
}
