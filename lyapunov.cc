#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <cmath>
#include <string>
#include <unistd.h> // for getopt()
#include <vector>
#include <omp.h>

#include "fractal_type.hh"
#include "string_case.hh"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

#define CLAMP(x) static_cast<uint8_t>(((x) > 255) ? (255u) : ((x) < 0) ? (0u) : (x))

//int seq_length = 12;
//int seq[12] = {1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0};

std::string center;
std::string seq = "BBBABBAAAAAA";
std::string type = "plain";
std::string format = "pnm";
std::string outfile;

double cmin = 2.5; // We go to 3rd dimension! \o/
double cmax = 4.0; //

double amin = 3.999;
double amax = 2.4;
double bmin = 2.4;
double bmax = 3.999;

double valueD = 3.7;
double valueE = 3.0;
double valueF = 3.0;

bool setC = false;
bool setD = false;
bool setE = false;
bool setF = false;

int asize = 800;
int bsize = 800;
int csize = 0;  // value >0 lets create 3df file instead of PPM
int gstart = 0; // first frame to render (0-based, of the csize frames)
int gend = 0;   // one past the last frame to render (0: render all frames)

int nmax = 1000;                /* number of rounds */

std::vector<double> params;

FILE* of = nullptr;
FractalType* fractal_type = nullptr;

/* for color generation; somewhat empirical, in order to match Wickerprint's original colors */
double lambda_min = -2.55;
double lambda_max = 0.3959;


std::vector<std::vector<float>> prob; // probability for each img pixel


struct RGB
{
	unsigned char r,g,b;
};


std::string replace_all(std::string orig, std::string_view search, std::string_view replacement)
{
	std::string::size_type pos = orig.find(search);
	while (pos < orig.size())
	{
		orig.replace(pos, search.length(), replacement);
		pos = orig.find(search, pos + replacement.length());
	}

	return orig;
}


std::string decorate_filename(std::string filename, const unsigned frame_nr)
{
	char frame_nr_string[16] = {};
	const unsigned digits = snprintf(nullptr, 0, "%u", csize-1);
	snprintf(frame_nr_string, sizeof(frame_nr_string), "%0*u", digits, frame_nr);
	filename = replace_all(filename, "{s}", seq);
	filename = replace_all(filename, "{i}", std::to_string(nmax) );
	filename = replace_all(filename, "{f}", frame_nr_string);
	return filename;
}

template<class PXL>
PXL lambda2pxl(double lambda);

template<>
RGB lambda2pxl<RGB>(const double lambda)
{
	double r=0,g=0,b=0;
	if (lambda > 0)
	{
		b = lambda/lambda_max;
	} else {
		r = 1 - pow(lambda/lambda_min, 2/3.0);
		g = 1 - pow(lambda/lambda_min, 1/3.0);
		b = 0;
	}

	return RGB{
		CLAMP(r*255),
		CLAMP(g*255),
		CLAMP(b*255)
	};
}

template<>
uint8_t lambda2pxl<uint8_t>(const double lambda)
{
	return CLAMP( 127 - lambda*100 );
}


// extension: 3rd dimension: value c is fixed, selected by 'C' in sequence string
template<class PXL>
void lyapunov(std::vector<PXL>& img, const int max_iter, const std::string& sequence, const FractalType& ftype, const unsigned frame)
{
	const int height = asize;
	const int width  = bsize;

	constexpr unsigned THREADS = 32;
	double tvals[THREADS][5] = { };

	for (unsigned t=0; t<THREADS; ++t)
	{
		ftype.per_image(tvals[t], frame);
	}

#pragma omp parallel for schedule(dynamic)
	for (int ai = 0; ai <height; ++ai)
	{
		const int tn = omp_get_thread_num();
		double* vals = tvals[tn];
//		const double a = amin + (amax-amin)/height*(ai+0.5);
//		vals[0] = amin + (amax-amin)/height*(ai+0.5);
		ftype.per_line(vals, frame, ai);

		fprintf(stderr, "\r%4d of %4d   ", ai, height);
		for (int bi = 0; bi < width; ++bi)
		{
//			const double b = bmin + (bmax-bmin)/width*(bi+0.5);
//			vals[1] = bmin + (bmax-bmin)/width*(bi+0.5);
			ftype.per_pixel(vals, frame, ai, bi);

			double x = 0.5;
			/* one round without derivating, so that the value 0.5 is avoided */
			for (unsigned m = 0; m < sequence.length(); m++)
			{
				const char ch = sequence[m];
				const double r = vals[ch-'A']; //ch=='B' ? b : (ch=='A' ? a : c);
				x = r*x*(1-x);
			}
			double sum_log_deriv = 0;
			double prod_deriv = 1;
			bool diverged = false;
			for (int n = 0; n < max_iter; n++)
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
				for (unsigned m = 0; m < sequence.length(); m++)
				{
					const char ch = sequence[m];
					const double r = vals[ch-'A']; // ch=='B' ? b : (ch=='A' ? a : c);
					round_deriv *= r*(1-2*x);
					x = r*x*(1-x);
				}
				/* accumulate the derivative product over several rounds and take
				   one logarithm only when it drifts out of a comfortable range
				   (one per round is what the old code did, one per pixel usually) */
				prod_deriv *= round_deriv;
				const double a = fabs(prod_deriv);
				if (a > 1e100 || a < 1e-100 || n == max_iter-1)
				{
					sum_log_deriv += log(a);
					prod_deriv = 1;
				}
			}
			const double lambda = diverged ? 1e30 : sum_log_deriv / (nmax*sequence.length());
			img[size_t(ai) * width + bi] = lambda2pxl<PXL>(lambda);
		}
	}
	
	fputs(" \n", stderr);
}


void make_ppm(const unsigned frame_nr)
{
	if (!outfile.empty())
	{
		std::string outfile_ppm = decorate_filename(outfile, frame_nr);
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
	img.resize(size_t(asize) * bsize);

	// PPM header
	fprintf( of,"P6\n");
	fflush(of);
	lyapunov(img, nmax, seq, *fractal_type, frame_nr);
	fprintf(of, "# Lyapunov: '%s' max_iter=%d  seq='%s'  a=%f ... %f  b=%f ... %f ",
		fractal_type->name(),
		nmax, seq.c_str(), amin, amax, bmin, bmax);

	if(setC) fprintf(of, " frame=%u/%u  C=%f ... %f", frame_nr, csize, cmin, cmax);
	if(setD) fprintf(of, " D=%f", valueD);
	if(setE) fprintf(of, " E=%f", valueE);
	if(setF) fprintf(of, " F=%f", valueF);
	
	fprintf(of, "\n"
		"%d %d 255\n",
		bsize, asize
		);
	
	for (int ai = 0; ai < asize; ++ai)
	{
		fwrite(&img[size_t(ai) * bsize], sizeof(RGB), bsize, of);
	}
	fflush(of);

	if (!outfile.empty())
	{
		fclose(of);
	}
}

void make_ppm()
{
	if (csize==0)
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
	if( unsigned(y) >= unsigned(asize) ) return 0;
	if( unsigned(x) >= unsigned(bsize) ) return 0;
	return img[z][size_t(y) * bsize + x];
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
		bsize >> 8, bsize & 255,
		asize >> 8, asize & 255,
		csize >> 8, csize & 255
	);
	fflush(of);
	
	// 3 images:
	std::vector<std::vector<uint8_t>> img;
	img.resize(3, std::vector<uint8_t>(size_t(asize) * bsize));
	
	for(int ci=0; ci<csize; ++ci)
	{
		const double c = (csize > 1) ? cmin + (cmax-cmin)/(csize-1)*ci : (cmin+cmax)/2;
		fprintf(stderr, "Frame %i of %i. c=%f \n", ci, csize, c);
		lyapunov(img[2], nmax, seq, *fractal_type, ci);
		
		fprintf(stderr,"\tblurring...\n");
		for(int y=0; y<asize; ++y)
		for(int x=0; x<bsize; ++x)
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
		"\t            the amplitude grows linearly from 0 to the given value over the frames (single frame: full amplitude)\n"
		"\tinversion : plain (A,B) plane, but the points are reflected about a fixed circle (P' = C + R^2*(P-C)/|P-C|^2)\n"
		"\t            -p [center_a,center_b[,radius]]  (default: window center, half of the A range)\n"
	);
}

void print_format_help(FILE* f)
{
	fprintf(f, "Output formats:\n"
		"\tppm : a PPM file or series of PPM files (if -o is given) or concatenated series of PPM images to stdout (of no -o given)\n"
		"\t3df : a 3DF voxel file\n"
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
int main(int argc, char** argv)
{
	const auto program_start = std::chrono::system_clock::now();
	int a;
	while( (a=getopt(argc, argv, "W:H:f:c:x:y:z:X:Y:Z:s:i:o:O:t:C:D:E:F:p:G:")) != -1 )
	{
		switch(a)
		{
			case 'W': bsize = atoi(optarg); break;
			case 'H': asize = atoi(optarg); break;
			case 'f': csize = atoi(optarg); setC = true; break;
			case 'z': cmin  = atof(optarg); setC = true; break;
			case 'Z': cmax  = atof(optarg); setC = true; break;
			case 'c': center = optarg;      break;
			case 'x': bmin = atof(optarg);  break;
			case 'y': amin = atof(optarg);  break;
			case 'X': bmax = atof(optarg);  break;
			case 'Y': amax = atof(optarg);  break;
			case 's': seq = optarg;         break;
			case 'i': nmax = atoi(optarg);  break;
			case 't': type = optarg;        break;
			case 'o': outfile = optarg;     break;
			case 'O': format = optarg;      break;
			case 'C': cmin   = atof(optarg); setC = true; break;
			case 'D': valueD = atof(optarg); setD = true; break;
			case 'E': valueE = atof(optarg); setE = true; break;
			case 'F': valueF = atof(optarg); setF = true; break;
			case 'p': params = parse_params(optarg); break;
			case 'G':
			{
				const char* colon = strchr(optarg, ':');
				if (colon == nullptr)
				{
					fprintf(stderr, "Option -G expects start:end.\n");
					return 2;
				}
				gstart = atoi(optarg);
				gend = atoi(colon+1);
				break;
			}
			case '?': [[fallthrough]];
			default:
				fprintf(stderr, 
					"Generate a Lyapunov fractal as PPM or 3df to stdout.\n\n"
					"Usage: %s [-W width] [-H height] [-f frames] [-s sequence] [-i max_iter]\n"
					"     [-t fractal_type ] [-o output_filename_base] [-O format]\n"
					"   { [-c center_x:center_y:size] | [-x min_x] [-X max_x] [-y min_y] [-Y max_y] } [-z min_z] [-Z max_z]\n"
					"     [-C value] [-D value] [-E value] [-F value]\n"
					"\n"
					"\t-W width  : width of image in pixel (default: %i)\n"
					"\t-H height : height of image in pixel (default: %i)\n"
					"\t-f frames : number of frames or size in in z direction\n"
					"\t-G start:end : render only the frames in the range [start, end) of the -f frames\n"
					"\t             (frame numbers and timing use the total -f count, useful for\n"
					"\t             rendering chunks in parallel; default: all frames)\n"
					"\t-s seq    : AB or ABC sequence for fractal (default: \"%s\")\n"
					"\t-i iter   : iterations per pixel (default: %i)\n"
					"\n"
					"\t-t type    : fractal type (default: %s, run \"-t help\" to get a list of types)\n"
					"\t-o outfile : basename of the output file (stdout if not given)\n"
					"\t             These patterns are replaced: {s}: sequence, {i}: iter, {f}: frame_nr\n"
					"\t-O format  : either \"pnm\" for image (series) or \"3df\" for 3D voxel file (default: \"%s\")\n"
					"\n"
					"\t-c cx:cy:x_size : center coordinates and size.\n"
					"\t   The x, X, y and Y values will be calculated from them automatically:\n"
					"\t-x min_x  : minimum x value (default: %f)\n"
					"\t-X max_x  : maximum x value (default: %f)\n"
					"\t-y min_y  : minimum y value (default: %f)\n"
					"\t-Y max_y  : maximum y value (default: %f)\n"
					"\t-z min_z  : 3rd dimension: (start) value of z. (selected by \"C\" in sequence string. default: %f)\n"
					"\t-Z max_z  : maximum value of z (default: %f) \n"
					"\t-C value  : value for 3th dimension (default %f)\n"
					"\t-D value  : value for 4th dimension (default %f)\n"
					"\t-E value  : value for 5th dimension (default %f)\n"
					"\n"
					"\t-p params[,...] : comma-separated list of additional numerical parameters\n"
					, argv[0], bsize, asize, seq.c_str(), nmax,
					type.c_str(), format.c_str(),
					bmin, bmax, amin, amax, cmin, cmax, cmin, valueD, valueE
					);
					return 1;
		}
	}
	
	if(csize > 0)
	{
		if(gend == 0) gend = csize; // no -G: render all frames
		if(gstart < 0 || gend <= gstart || gend > csize)
		{
			fprintf(stderr, "Invalid frame range -G %d:%d (total frames: %d).\n", gstart, gend, csize);
			return 2;
		}
	}

	if(seq.empty())
	{
		fprintf(stderr, "Sequence must not be empty.\n");
		return 2;
	}

	const auto illegal = std::find_if(seq.begin(), seq.end(), [](char ch) { return ch < 'A' || ch > 'E'; });
	if (illegal != seq.end())
	{
		fprintf(stderr, "Sequence \"%s\" contains illegal character '%c' (only 'A'..'E' are allowed).\n", seq.c_str(), *illegal);
		return 2;
	}

	if(nmax <= 0)
	{
		fprintf(stderr, "max_iter must be > 0 (got %d).\n", nmax);
		return 2;
	}

	if(!center.empty())
	{
		double cx=0.0, cy=0.0, dx=0.0;
		const int ret = sscanf(center.c_str(), "%lg:%lg:%lg", &cx, &cy, &dx);
		if(ret != 3)
		{
			fprintf(stderr, "Cannot parse \"%s\" as \"cx:cy:x_size\" parameter.\n", center.c_str());
			return 2;
		}
		
		const double dy = (dx * asize) / bsize;
		bmin = cx - dx/2;
		bmax = cx + dx/2;
		amin = cy + dy/2;
		amax = cy - dy/2;
		
		fprintf(stderr, "Interval calculated from \"-c %s\" ", center.c_str());
	}

	fprintf(stderr, " min_x=%11.9f  max_x=%11.9f  min_y=%11.9f  max_y=%11.9f", bmin, bmax, amin, amax);

	if(setC) fprintf(stderr, " C=%f ... %f", cmin, cmax);
	if(setD) fprintf(stderr, " D=%f", valueD);
	if(setE) fprintf(stderr, " E=%f", valueE);
	if(setF) fprintf(stderr, " F=%f", valueF);
	if (!params.empty())
	{
		fprintf(stderr, " P=%f", params[0]);
		for (unsigned u=1; u<params.size(); ++u)
		{
			fprintf(stderr, ", %f", params[u]);
		}
	}

	fputc('\n', stderr);


	fprintf(stderr, "Run with %d threads...\n", omp_get_max_threads());
	of = stdout;

#ifdef _WIN32
	setmode(fileno(stdout),O_BINARY);
	setmode(fileno(stdin),O_BINARY);
#endif

	const Dimension dim_a{unsigned(asize), amin, amax};
	const Dimension dim_b{unsigned(bsize), bmin, bmax};
	const Dimension dim_c{unsigned(csize), cmin, cmax};

	switch (case_hash(type))
	{
		case "plain"_case  : fractal_type = new SimpleFractalType(dim_a, dim_b, dim_c, valueD, valueE); break;
		case "6sides"_case : fractal_type = new Flight(dim_a, dim_b, dim_c, valueD, valueE, std::make_unique<SixSidePath>(params)); break;
		case "wander"_case : fractal_type = new Flight(dim_a, dim_b, dim_c, valueD, valueE, std::make_unique<WanderPath>(params)); break;
		case "circle"_case : fractal_type = new Circle(dim_a, dim_b, dim_c, valueD, valueE); break;
		case "sphere"_case : fractal_type = new Sphere(dim_a, dim_b, dim_c, valueD, valueE, params); break;
		case "tileable"_case : fractal_type = new Tileable(dim_a, dim_b, dim_c, valueD, valueE); break;
		case "tileable-loop"_case : fractal_type = new TileableLoop(dim_a, dim_b, dim_c, valueD, valueE); break;
		case "zoom"_case : fractal_type = new Zoom(dim_a, dim_b, dim_c, valueD, valueE, params); break;
		case "wavy"_case :
			if (params.size() > 0 && params[0] <= 0.0)
			{
				fprintf(stderr, "wavy: frequency_a must be > 0 (got %f).\n", params[0]);
				return 2;
			}
			if (params.size() > 3 && params[3] <= 0.0)
			{
				fprintf(stderr, "wavy: frequency_b must be > 0 (got %f).\n", params[3]);
				return 2;
			}
			fractal_type = new Wavy(dim_a, dim_b, dim_c, valueD, valueE, params); break;
		case "noise"_case :
			if (params.size() > 0 && params[0] <= 0.0)
			{
				fprintf(stderr, "noise: frequency must be > 0 (got %f).\n", params[0]);
				return 2;
			}
			fractal_type = new Noise(dim_a, dim_b, dim_c, valueD, valueE, params); break;
		case "inversion"_case : fractal_type = new Inversion(dim_a, dim_b, dim_c, valueD, valueE, params); break;
		
		case "help"_case   : [[fallthrough]];
		case "-h"_case     : [[fallthrough]];
		case "--help"_case : print_type_help(stdout); return 0;
		default: print_type_help(stderr); return 1;
	}

	if (case_hash(format) == "3df"_case && csize <= 0)
	{
		fprintf(stderr, "3df output requires a z dimension: use -f <frames>.\n");
		return 2;
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
