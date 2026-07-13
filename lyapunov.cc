#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
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

#define CLAMP(x) (((x) > 255) ? (unsigned char)255u : ((x) < 0) ? (unsigned char)0u : (unsigned char)(x))

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

int nmax = 1000;                /* number of rounds */

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
void lyapunov(std::vector<std::vector<PXL>>& img, const int max_iter, const std::string& sequence, const FractalType& ftype, const unsigned frame)
{
	double min_x = +HUGE_VAL;
	double max_x = -HUGE_VAL;
	double min_lambda = +HUGE_VAL;
	double max_lambda = -HUGE_VAL;

	const int height = img.size();
	const int width  = img[0].size();

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
			for (int n = 0; n < max_iter; n++)
			{
				double prod_deriv = 1;
				for (unsigned m = 0; m < sequence.length(); m++)
				{
					const char ch = sequence[m];
					const double r = vals[ch-'A']; // ch=='B' ? b : (ch=='A' ? a : c);
					/* avoid computing too many logarithms. One every round is acceptable. */
					prod_deriv *= r*(1-2*x); 
					x = r*x*(1-x);
				}
				sum_log_deriv += log(fabs(prod_deriv));
			}
			const double lambda = sum_log_deriv / (nmax*sequence.length());
			/*
			if(x>max_x) max_x = x;
			if(x<min_x) min_x = x;
			if(lambda>max_lambda) max_lambda = lambda;
			if(lambda<min_lambda) min_lambda = lambda;
			*/
			img[ai][bi] = lambda2pxl<PXL>(lambda);
		}
	}
	
	fprintf(stderr, " x=%f ... %f ; lambda=%f ... %f \n",
			min_x, max_x, min_lambda, max_lambda );
	
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
	std::vector<std::vector<RGB>> img;
	img.resize(asize);
	for(auto& row : img)
	{
		row.resize(bsize);
	}

	// PPM header
	fprintf( of,"P6\n");
	fflush(of);
	lyapunov(img, nmax, seq, *fractal_type, frame_nr);
	fprintf(of, "# Lyapunov: '%s' max_iter=%d  seq='%s'  a=%f ... %f  b=%f ... %f ",
		fractal_type->name(),
		nmax, seq.c_str(), amin, amax, bmin, bmax);
	
	if(setD) fprintf(of," D=%f", valueD);
	if(setE) fprintf(of, " E=%f", valueE);
	if(setF) fprintf(of, " F=%f", valueF);
	
	fprintf(of, "\n"
		"%d %d 255\n",
		bsize, asize
		);
	
	for(const auto& row : img)
	{
		fwrite(row.data(), sizeof(RGB), row.size(), of);
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
	}else for (int frame_nr=0; frame_nr < csize; ++frame_nr)
	{
		make_ppm(frame_nr);
	}
}

inline
uint8_t get(const std::vector<std::vector<std::vector<uint8_t>>>& img, const int x, const int y, const int z)
{
//	if( unsigned(z)>=img.size() ) return 0;
	if( unsigned(y)>=img[z].size() ) return 0;
	if( unsigned(x)>=img[z][y].size() ) return 0;
	return img[z][y][x];
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
	std::vector<std::vector<std::vector<uint8_t>>> img;
	img.resize(3);
	img[0].resize(asize);
	for(auto& row : img[0])
	{
		row.resize(bsize);
	}
	
	img[2] = img[1] = img[0];
	
	for(int ci=0; ci<csize; ++ci)
	{
		const double c = cmin + (cmax-cmin)/(csize-1)*ci;
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
				u += get(img, x+dx, y+dy, z );
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
		"\t6sides : animated cruise through all 6 sides of the ABC cube. D=fix, E=fix\n"
		"\tcircle : C+D circle\n"
	);
}

void print_format_help(FILE* f)
{
	fprintf(f, "Output formats:\n"
		"\tppm : a PPM file or series of PPM files (if -o is given) or concatenated series of PPM images to stdout (of no -o given)\n"
		"\t3df : a 3DF voxel file\n"
	);
}

int main(int argc, char** argv)
{
	const auto program_start = std::chrono::system_clock::now();
	int a;
	while( (a=getopt(argc, argv, "W:H:f:c:x:y:z:X:Y:Z:s:i:o:O:t:C:D:E:F:")) != -1 )
	{
		switch(a)
		{
			case 'W': bsize = atoi(optarg); break;
			case 'H': asize = atoi(optarg); break;
			case 'f': csize = atoi(optarg); break;
			case 'z': cmin  = atof(optarg); break;
			case 'Z': cmax  = atof(optarg); break;
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
					"\t-s seq    : AB or ABC sequence for fractal (default: \"%s\")\n"
					"\t-i iter   : iterations per pixel (default: %i)\n"
					"\n"
					"\t-t type    : fractal type (default: %s, run \"-t help\" to get a list of types)\n"
					"\t-o outfile : basename of the output file (stdout if not given)\n"
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
					, argv[0], bsize, asize, seq.c_str(), nmax,
					type.c_str(), format.c_str(),
					bmin, bmax, amin, amax, cmin, cmax, cmin, valueD, valueE
					);
					return 1;
		}
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
		
		fprintf(stderr, "Interval calculated from \"-c %s\": min_x=%11.9f  max_x=%11.9f  min_y=%11.9f  max_y=%11.9f",
			center.c_str(), bmin, bmax, amin, amax);
		
		if(setC) fprintf(stderr, " C=%f", cmin);
		if(setD) fprintf(stderr, " D=%f", valueD);
		if(setE) fprintf(stderr, " E=%f", valueE);
		if(setF) fprintf(stderr, " F=%f", valueF);
		
		fputc('\n', stderr);
	}

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
		case "6sides"_case : fractal_type = new SixSide(dim_a, dim_b, dim_c, valueD, valueE); break;
		case "circle"_case : throw std::logic_error("Unimplemented"); break;
		
		case "help"_case   : [[fallthrough]];
		case "-h"_case     : [[fallthrough]];
		case "--help"_case : print_type_help(stdout); return 0;
		default: print_type_help(stderr); return 1;
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
