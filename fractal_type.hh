
#ifndef FRACTAL_TYPE_HH
#define FRACTAL_TYPE_HH

#include <array>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Dimension
{
	unsigned size; // in pixel
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
	struct Quat { double w, x, y, z; };
	struct Vec3 { double x, y, z; };

	// the (A,B) slice at C=cmin tilts about 6 cube edges (90 degrees each,
	// axes A,C,B,A,C,B), visiting all 6 faces of the ABC cube
	// and ending back at the start face (window rotated 180 degrees)
	// _params: [rotate][,ease]
	//   rotate: 1: the window additionally rotates 180 degrees in the plane
	//           over the whole animation, so the animation loops seamlessly
	//   ease:   1: each 90-degree tilt uses smoothstep timing (3t^2-2t^3)
	//           instead of a constant angular speed, so the window settles
	//           at every face and the rotation axis switches are not visible
	//           2: quaternion Catmull-Rom path through the 6 face poses: the
	//           rotation axis blends continuously between the tilts (cruise,
	//           no pause at the faces)
	//   default: 0,0
	SixSide(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, rotate(_params.size() > 0 && _params[0] != 0.0)
	, ease(_params.size() > 1 ? (int)_params[1] : 0)
	{
		if (ease == 2) build_cr_path();
	}

	void per_image(double values[5], unsigned frame) const override;
	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "six_side"; }

private:
	const bool rotate;
	const int ease;
	// ease == 2: periodic Catmull-Rom path through the 7 face poses K0..K6
	// (K6 = start face rotated 180 degrees); cr_qend[0] is the pose at t=-1
	// (= K5), cr_qend[1] the pose at t=7 (= K1), with aligned quaternion signs
	std::array<Quat, 7> cr_qk{};
	std::array<Vec3, 7> cr_tk{};
	std::array<Quat, 2> cr_qend{};
	std::array<Vec3, 2> cr_tend{};
	// cumulative path length for the constant-speed parameterization
	std::vector<double> cr_arclen;
	void build_cr_path();
	void cr_pose_at(double t, Quat& q, Vec3& tr) const;
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

class Wavy : public SimpleFractalType
{
public:
	// _params: frequency_a (waves per B axis), amplitude_a (ratio of image size), phase_a (degrees),
	//          [frequency_b, amplitude_b, phase_b] (B defaults to the A values)
	//          default triple: 1, 0.02, 0 (360 degrees = one full wave)
	Wavy(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, freq_a(_params.size()>0 ? _params[0] : 1.0)
	, amp_a_ratio(_params.size()>1 ? _params[1] : 0.02)
	, phase_a(2.0*M_PI * (_params.size()>2 ? _params[2] : 0.0) / 360.0)
	, freq_b(_params.size()>3 ? _params[3] : freq_a)
	, amp_b_ratio(_params.size()>4 ? _params[4] : amp_a_ratio)
	, phase_b(_params.size()>5 ? 2.0*M_PI*_params[5]/360.0 : phase_a)
	, amp_a(amp_a_ratio * (a.max-a.min))
	, amp_b(amp_b_ratio * (b.max-b.min))
	, k_a(2.0*M_PI*freq_a/b.size)
	, k_b(2.0*M_PI*freq_b/a.size)
	{}

	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "wavy"; }

private:
	const double freq_a, amp_a_ratio, phase_a;
	const double freq_b, amp_b_ratio, phase_b;
	// phases in radians (input in degrees)
	// amplitude as fraction of the A/B parameter range
	const double amp_a, amp_b;
	// angular frequency, radians per pixel (A wobbles along columns, B along lines)
	const double k_a, k_b;
};

namespace detail
{
// 2D Perlin gradient noise, value roughly in [-1, 1]
struct Perlin
{
	explicit Perlin(unsigned seed);
	[[nodiscard]] double noise(double x, double y) const;

private:
	int perm[512] = {};
};
}

class Noise : public SimpleFractalType
{
public:
	// _params: frequency (noise cells per image axis), amplitude_a, amplitude_b (ratio of image size),
	//          [seed] (default 3, 0.333, 0.333, 0)
	// the amplitude grows linearly from 0 to the given value over the frames (single frame: full amplitude)
	Noise(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, freq(_params.size()>0 ? _params[0] : 3.0)
	, amp_a((_params.size()>1 ? _params[1] : 0.333) * (a.max-a.min))
	, amp_b((_params.size()>2 ? _params[2] : 0.333) * (b.max-b.min))
	, perlin(_params.size()>3 ? (unsigned)_params[3] : 0)
	{}

	void per_line (double values[5], unsigned frame, unsigned line) const override;
	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "noise"; }

private:
	const double freq, amp_a, amp_b;
	const detail::Perlin perlin;
};

class Inversion : public SimpleFractalType
{
public:
	// _params: [center_a, center_b] (default window center),
	//          [radius] (default half of the A range, i.e. "half the image size")
	Inversion(const Dimension& _a, const Dimension& _b, const Dimension& _c, double _d, double _e, const std::vector<double>& _params)
	: SimpleFractalType(_a, _b, _c, _d, _e)
	, center_a(_params.size()>0 ? _params[0] : (a.max+a.min)/2)
	, center_b(_params.size()>1 ? _params[1] : (b.max+b.min)/2)
	, radius(_params.size()>2 ? _params[2] : (a.max-a.min)/2)
	{}

	void per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const override;
	[[nodiscard]] const char* name() const override { return "inversion"; }

private:
	const double center_a, center_b, radius;
};
#endif // FRACTAL_TYPE_HH
