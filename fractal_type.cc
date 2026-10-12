#include "fractal_type.hh"

#include <math.h>
#include <algorithm>

namespace
{
using Quat = SixSide::Quat;
using Vec3 = SixSide::Vec3;

Quat qscale(Quat a, double s) { return {a.w*s, a.x*s, a.y*s, a.z*s}; }
Quat qadd(Quat a, Quat b) { return {a.w+b.w, a.x+b.x, a.y+b.y, a.z+b.z}; }
double qdot(Quat a, Quat b) { return a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z; }
Quat qnormu(Quat a) { return qscale(a, 1.0/std::sqrt(qdot(a, a))); }

Vec3 cross3(Vec3 a, Vec3 b)
{
	return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

// rotation matrix (columns: u-axis, v-axis, normal) to unit quaternion
Quat mat2q(const Vec3& ru, const Vec3& rv, const Vec3& n)
{
	const double m00 = ru.x, m01 = ru.y, m02 = ru.z;
	const double m10 = rv.x, m11 = rv.y, m12 = rv.z;
	const double m20 = n.x,  m21 = n.y,  m22 = n.z;
	const double tr = m00 + m11 + m22;
	Quat q;
	if (tr > 0.0)
	{
		const double S = std::sqrt(tr + 1.0) * 2.0;
		q.w = 0.25*S;    q.x = (m12 - m21)/S; q.y = (m20 - m02)/S; q.z = (m01 - m10)/S;
	}
	else if (m00 > m11 && m00 > m22)
	{
		const double S = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
		q.w = (m12 - m21)/S; q.x = 0.25*S;    q.y = (m01 + m10)/S; q.z = (m20 + m02)/S;
	}
	else if (m11 > m22)
	{
		const double S = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
		q.w = (m20 - m02)/S; q.x = (m01 + m10)/S; q.y = 0.25*S;    q.z = (m21 + m12)/S;
	}
	else
	{
		const double S = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
		q.w = (m01 - m10)/S; q.x = (m20 + m02)/S; q.y = (m21 + m12)/S; q.z = 0.25*S;
	}
	return qnormu(q);
}

Vec3 qrot(Quat q, Vec3 v)
{
	const Vec3 t{2.0*(q.y*v.z - q.z*v.y), 2.0*(q.z*v.x - q.x*v.z), 2.0*(q.x*v.y - q.y*v.x)};
	return {v.x + q.w*t.x + (q.y*t.z - q.z*t.y),
	        v.y + q.w*t.y + (q.z*t.x - q.x*t.z),
	        v.z + q.w*t.z + (q.x*t.y - q.y*t.x)};
}

// Catmull-Rom (Hermite) basis on [0,1]
double cr_h00(double t) { return (2.0*t - 3.0)*t*t + 1.0; }     // 2t^3-3t^2+1
double cr_h01(double t) { return t*t*(3.0 - 2.0*t); }           // -2t^3+3t^2
double cr_ht(double t)  { return t*t*(2.0*t - 2.0); }           // (h10+h11) = 2t^3-2t^2

Quat cr_quat(const Quat& qm, const Quat& q0, const Quat& q1, const Quat& qp, double t)
{
	return qnormu(qadd(qadd(qscale(q0, cr_h00(t)), qscale(q1, cr_h01(t))),
	                   qscale(qadd(qp, qscale(qm, -1.0)), 0.5*cr_ht(t))));
}

Vec3 cr_vec3(const Vec3& tm, const Vec3& t0, const Vec3& t1, const Vec3& tp, double t)
{
	const double a = cr_h00(t), b = cr_h01(t), c = 0.5*cr_ht(t);
	const Vec3 d{tp.x - tm.x, tp.y - tm.y, tp.z - tm.z};
	return {t0.x*a + t1.x*b + d.x*c, t0.y*a + t1.y*b + d.y*c, t0.z*a + t1.z*b + d.z*c};
}

// the 6 hinge segments: (A01,B01,C01) from window coordinates (u,v) at tilt theta
Vec3 sixside_map(unsigned seg, double theta, double u, double v)
{
	const double s = std::sin(theta), co = std::cos(theta);
	switch (seg)
	{
	case 0:  return {v, 1-(1-u)*co, (1-u)*s};
	case 1:  return {1-(1-v)*co, 1-(1-v)*s, 1-u};
	case 2:  return {1-u*s, v, 1-u*co};
	case 3:  return {1-u, v*co, 1-v*s};
	case 4:  return {(1-u)*co, (1-u)*s, 1-v};
	default: return {(1-v)*s, 1-u, (1-v)*co};
	}
}
} // namespace

void SixSide::build_cr_path()
{
	for (unsigned k = 0; k < 7; ++k)
	{
		// K0..K5: segment k at tilt 0, K6: segment 5 at tilt 90 degrees
		const unsigned seg = (k < 6) ? k : 5;
		const double theta = (k == 6) ? M_PI/2.0 : 0.0;
		const Vec3 p0 = sixside_map(seg, theta, 0.0, 0.0);
		const Vec3 pu = sixside_map(seg, theta, 1.0, 0.0);
		const Vec3 pv = sixside_map(seg, theta, 0.0, 1.0);
		const Vec3 ru{pu.x-p0.x, pu.y-p0.y, pu.z-p0.z};
		const Vec3 rv{pv.x-p0.x, pv.y-p0.y, pv.z-p0.z};
		cr_tk[k] = p0;
		cr_qk[k] = mat2q(ru, rv, cross3(ru, rv));
	}
	for (unsigned k = 1; k < 7; ++k)
		if (qdot(cr_qk[k-1], cr_qk[k]) < 0.0)
			cr_qk[k] = qscale(cr_qk[k], -1.0);
	cr_qend[0] = qdot(cr_qk[5], cr_qk[0]) < 0.0 ? qscale(cr_qk[5], -1.0) : cr_qk[5];
	cr_qend[1] = qdot(cr_qk[1], cr_qk[6]) < 0.0 ? qscale(cr_qk[1], -1.0) : cr_qk[1];
	cr_tend[0] = cr_tk[5];
	cr_tend[1] = cr_tk[1];

	// cumulative path length over the whole tour, used to reparameterize the
	// path by constant speed (the raw Catmull-Rom timing accelerates and
	// decelerates between the faces, which looks like a stutter)
	const int m = 6 * 512;
	cr_arclen.resize(m + 1);
	cr_arclen[0] = 0.0;
	auto [q0, t0] = cr_pose_at(0.0);
	for (int i = 1; i <= m; ++i)
	{
		const auto [q, tr] = cr_pose_at(6.0 * i / m);
		const double d = std::min(1.0, std::abs(qdot(q0, q)));
		const double ang = 2.0 * std::acos(d);
		const double dt = std::sqrt((tr.x-t0.x)*(tr.x-t0.x) + (tr.y-t0.y)*(tr.y-t0.y) + (tr.z-t0.z)*(tr.z-t0.z));
		cr_arclen[i] = cr_arclen[i-1] + std::sqrt(ang*ang + dt*dt);
		q0 = q;
		t0 = tr;
	}
}

std::pair<Quat, Vec3> SixSide::cr_pose_at(double t) const
{
	const unsigned seg = std::min(t, 5.9999999);
	const double tau = t - seg;
	const Quat qm = (seg == 0) ? cr_qend[0] : cr_qk[seg-1];
	const Quat qp = (seg == 5) ? cr_qend[1] : cr_qk[seg+2];
	const Vec3 tm = (seg == 0) ? cr_tend[0] : cr_tk[seg-1];
	const Vec3 tp = (seg == 5) ? cr_tend[1] : cr_tk[seg+2];
	return {cr_quat(qm, cr_qk[seg], cr_qk[seg+1], qp, tau),
	        cr_vec3(tm, cr_tk[seg], cr_tk[seg+1], tp, tau)};
}

void SixSide::per_image(double values[5], unsigned frame) const
{
	values[3] = d;
	values[4] = e;
}

void SixSide::per_line(double values[5], unsigned frame, unsigned line) const
{
	// A depends on line and column, see per_pixel
}

void SixSide::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	double u = (column+0.5) / b.size;
	double v = (line+0.5) / a.size;
	const double t = (c.size > 1) ? 6.0 * (frame+0.5) / c.size : 0.0;
	const unsigned seg = std::min(t, 5.9999999);
	const double tau = t - seg;

	if (rotate && c.size > 1)
	{
		// additional in-plane rotation of the window, 0..180 degrees over
		// the whole animation, so that the end pose equals the start pose
		const double psi = (frame + 0.5) * M_PI / c.size;
		const double ps = std::sin(psi), pc = std::cos(psi);
		const double du = u - 0.5, dv = v - 0.5;
		u = du*pc - dv*ps + 0.5;
		v = du*ps + dv*pc + 0.5;
	}

	// normalized cube coordinates (A01,B01,C01) in [0,1]^3, material window
	// coordinates u = column ratio, v = line ratio
	double A01, B01, C01;
	if (ease == 2)
	{
		// quaternion Catmull-Rom path through the 7 face poses (periodic,
		// so the path is also smooth across the animation loop point),
		// reparameterized by arc length for a constant cruise speed
		const int n = cr_arclen.size() - 1;
		const double target = (t / 6.0) * cr_arclen[n];
		int lo = 0, hi = n;
		while (hi - lo > 1)
		{
			const int mid = (lo + hi) / 2;
			if (cr_arclen[mid] < target) lo = mid; else hi = mid;
		}
		const double span = cr_arclen[hi] - cr_arclen[lo];
		const double frac = (span > 0.0) ? (target - cr_arclen[lo]) / span : 0.0;
		const double path_t = (lo + frac) * 6.0 / n;
		const auto [q, tr] = cr_pose_at(path_t);
		const Vec3 p = qrot(q, {u, v, 0.0});
		A01 = p.x + tr.x;  B01 = p.y + tr.y;  C01 = p.z + tr.z;
	}
	else
	{
		const double theta = (M_PI/2.0) * (ease ? tau*tau*(3.0 - 2.0*tau) : tau);
		const double s = std::sin(theta), co = std::cos(theta);
		switch (seg)
		{
		case 0:  // face C=cmin -> B=bmax, hinged at the (B=bmax,C=cmin) edge, axis A
			A01 = v;              B01 = 1-(1-u)*co;   C01 = (1-u)*s;         break;
		case 1:  // face B=bmax -> A=amax, hinged at the (A=amax,B=bmax) edge, axis C
			A01 = 1-(1-v)*co;     B01 = 1-(1-v)*s;    C01 = 1-u;             break;
		case 2:  // face A=amax -> C=cmax, hinged at the (A=amax,C=cmax) edge, axis B
			A01 = 1-u*s;          B01 = v;            C01 = 1-u*co;          break;
		case 3:  // face C=cmax -> B=bmin, hinged at the (B=bmin,C=cmax) edge, axis A
			A01 = 1-u;            B01 = v*co;         C01 = 1-v*s;           break;
		case 4:  // face B=bmin -> A=amin, hinged at the (A=amin,B=bmin) edge, axis C
			A01 = (1-u)*co;       B01 = (1-u)*s;      C01 = 1-v;             break;
		default: // face A=amin -> C=cmin (start face, window rotated 180 degrees), hinged at the (A=amin,C=cmin) edge, axis B
			A01 = (1-v)*s;        B01 = 1-u;          C01 = (1-v)*co;        break;
		}
	}

	values[0] = a.min + A01 * (a.max-a.min);
	values[1] = b.min + B01 * (b.max-b.min);
	values[2] = c.min + C01 * (c.max-c.min);
}

void Circle::per_image(double values[5], unsigned frame) const
{
	if (c.size==0)
	{
		values[2] = c.min;
		values[3] = d;
		values[4] = e;
	}else{
		const double cm = (c.min + c.max)/2;
		const double r = abs(c.max - c.min);
		const double phi = frame * 2*M_PI / c.size;
		values[2] = cm - r*cos(phi);
		values[3] = d + r*sin(phi);
		values[4] = e;
	}
}


void Sphere::per_image(double values[5], unsigned frame) const
{
	const double phi = (c.size > 0) ? frame * 2*M_PI / c.size : 0.0;
	values[3] = d - radius_de * cos(phi);
	values[4] = e + radius_de * sin(phi);
}

void Sphere::per_line(double values[5], unsigned frame, unsigned line) const
{
	values[0] = center_a - sc_line.at(line).c * radius_a;
}

void Sphere::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const auto pl = sc_line.at(line);
	const auto pc = sc_column.at(column);
	values[1] = center_b + pl.s * pc.c * radius_b;
	values[2] = center_b - pl.s * pc.s * radius_c;
}

std::vector<SinCos> Sphere::calculate_line_vector() const
{
	std::vector<SinCos> result;
	result.reserve(a.size);
	for (unsigned i = 0; i < a.size; ++i)
	{
		const double phi = i*M_PI/a.size;
		const double ss = sin(phi);
		const double cc = cos(phi);
		result.emplace_back(ss,cc);
	}
	return result;
}

std::vector<SinCos> Sphere::calculate_column_vector() const
{
	std::vector<SinCos> result;
	result.reserve(b.size);
	for (unsigned i = 0; i < b.size; ++i)
	{
		const double phi = i*(2*M_PI)/b.size;
		const double ss = sin(phi);
		const double cc = cos(phi);
		result.emplace_back(ss,cc);
	}
	return result;
}

void Tileable::per_line(double values[5], unsigned frame, unsigned line) const
{
	values[0] = center_a + radius_a * sin_line.at(line);
}

void Tileable::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	values[1] = center_b + radius_b * cos_column.at(column);
}

std::vector<double> Tileable::calculate_sin_line() const
{
	std::vector<double> result;
	result.reserve(a.size);
	for (unsigned i = 0; i < a.size; ++i)
		result.push_back(sin(i*(2*M_PI)/a.size));
	return result;
}

std::vector<double> Tileable::calculate_cos_column() const
{
	std::vector<double> result;
	result.reserve(b.size);
	for (unsigned i = 0; i < b.size; ++i)
		result.push_back(cos(i*(2*M_PI)/b.size));
	return result;
}

void TileableLoop::per_image(double values[5], unsigned frame) const
{
	if (c.size==0)
	{
		values[2] = c.min;
		values[3] = d;
		values[4] = e;
	}else{
		const double cm = (c.min + c.max)/2;
		const double r = abs(c.max - c.min);
		const double phi = frame * 2*M_PI / c.size;
		values[2] = cm - r*cos(phi);
		values[3] = d + r*sin(phi);
		values[4] = e;
	}
}

void Zoom::per_image(double values[5], unsigned frame) const
{
	if (frame != cache_frame)
	{
		const double span = (c.size > 1) ? (double)(c.size - 1) : 1.0;
		const double t = (double)frame / span;
		const double s = pow(zoom_total, t);
		const double theta = rot_total * t;
		a_x = (a.max-a.min) * s * sin(theta);
		a_y = (a.max-a.min) * s * cos(theta);
		b_x = (b.max-b.min) * s * cos(theta);
		b_y = (b.max-b.min) * s * sin(theta);
		cache_frame = frame;
	}
	values[2] = center_c;
	values[3] = d;
	values[4] = e;
}

void Zoom::per_line(double values[5], unsigned frame, unsigned line) const
{
}

void Zoom::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const double x = (column + 0.5) / b.size - 0.5;
	const double y = (line + 0.5) / a.size - 0.5;
	values[0] = center_a + a_x * x + a_y * y;
	values[1] = center_b + b_x * x - b_y * y;
}

void Wavy::per_line(double values[5], unsigned frame, unsigned line) const
{

}

void Wavy::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5)   + amp_a * sin(k_a * column - phase_a);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5) + amp_b * sin(k_b * line - phase_b);
}

detail::Perlin::Perlin(unsigned seed)
{
	int p[256];
	for (int i = 0; i < 256; ++i)
		p[i] = i;
	unsigned int s = seed ? seed : 1;
	for (int i = 255; i > 0; --i)
	{
		s = s*1103515245u + 12345u;
		const int j = (int)((s >> 16) % (unsigned)(i+1));
		std::swap(p[i], p[j]);
	}
	for (int i = 0; i < 512; ++i)
		perm[i] = p[i & 255];
}

double detail::Perlin::noise(double x, double y) const
{
	const int xi = (int)floor(x) & 255;
	const int yi = (int)floor(y) & 255;
	const double xf = x - floor(x);
	const double yf = y - floor(y);
	const double u = xf*xf*xf*(xf*(xf*6-15)+10);
	const double v = yf*yf*yf*(yf*(yf*6-15)+10);

	auto dot2 = [](int hash, double dx, double dy) -> double
	{
		switch (hash & 7)
		{
		case 0: return  dx + dy;
		case 1: return  dx - dy;
		case 2: return -dx + dy;
		case 3: return -dx - dy;
		case 4: return  dx;
		case 5: return -dx;
		case 6: return  dy;
		default: return -dy;
		}
	};

	const double aa = dot2(perm[perm[xi] + yi],     xf,     yf);
	const double ba = dot2(perm[perm[xi+1] + yi],   xf-1,   yf);
	const double ab = dot2(perm[perm[xi] + yi+1],   xf,     yf-1);
	const double bb = dot2(perm[perm[xi+1] + yi+1], xf-1,   yf-1);
	const double x1 = aa + u*(ba - aa);
	const double x2 = ab + u*(bb - ab);
	return x1 + v*(x2 - x1);
}

void Noise::per_line(double values[5], unsigned frame, unsigned line) const
{
	// A depends on both line and column (2D noise)
}

void Noise::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	const double x = (column + 0.5) / b.size * freq;
	const double y = (line + 0.5) / a.size * freq;
	const double t = (c.size > 1) ? (double)frame / (c.size - 1) : 1.0;
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5) + amp_a * t * perlin.noise(x, y);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5) + amp_b * t * perlin.noise(x, y + 52.0);
}

void Inversion::per_pixel(double values[5], unsigned frame, unsigned line, unsigned column) const
{
	// reflection about the inversion circle: distance d and angle theta from the center,
	// projected to the radius R^2/d (radius R scaled by R/d) at the same angle
	values[0] = a.min + (a.max-a.min)/a.size*(line+0.5);
	values[1] = b.min + (b.max-b.min)/b.size*(column+0.5);
	const double da = values[0] - center_a;
	const double db = values[1] - center_b;
	const double d = sqrt(da*da + db*db);
	if (d > 1e-150) // the center itself is fixed (its image is at infinity)
	{
		const double scale = radius*radius / (d*d);
		values[0] = center_a + scale * da;
		values[1] = center_b + scale * db;
	}
}


