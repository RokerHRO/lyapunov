#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <memory>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>
#include <getopt.h>
#include <utility>
#include <vector>
#include <omp.h>

#include "fractal_type.hh"
#include "string_case.hh"
#include "rgb.hh"
#include "draw.hh"
#include "grid.hh"
#include "file_writer.hh"
#include "render.hh"

#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif

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
GridOptions grid_options; // set in main() from grid_spec

FILE* of = nullptr;
std::unique_ptr<FractalType> fractal_type;



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

	render_lyapunov(img, fractal_params, *fractal_type, frame_nr);
	if (grid)
		draw_grid(img, fractal_params, *fractal_type, !warp_name.empty(), grid_options);

	// PPM comment with the generation parameters
	std::ostringstream oss;
	oss << std::fixed << std::setprecision(6);
	oss << "# Lyapunov: '" << fractal_type->name() << "' max_iter=" << fractal_params.max_iter
	    << "  seq='" << fractal_params.sequence << "'  a=" << fractal_params.a.min << " ... " << fractal_params.a.max
	    << "  b=" << fractal_params.b.min << " ... " << fractal_params.b.max;
	// the optional parts below start with a space; the base always ends with
	// one (as in the original format string), so two spaces separate them
	oss << " ";
	if(setC) oss << " frame=" << frame_nr << "/" << fractal_params.c.size << "  C=" << fractal_params.c.min << " ... " << fractal_params.c.max;
	if(setD) oss << " D=" << fractal_params.d;
	if(setE) oss << " E=" << fractal_params.e;
	if(setF) oss << " F=" << valueF;

	write_ppm(of, img, fractal_params.b.size, fractal_params.a.size, oss.str());

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
	write_3df_header(of, fractal_params.b.size, fractal_params.a.size, fractal_params.c.size);
	fflush(of);
	
	// 3 images:
	std::vector<std::vector<uint8_t>> img;
	img.resize(3, std::vector<uint8_t>(fractal_params.a.size * fractal_params.b.size));
	
	for(int ci=0; ci<int(fractal_params.c.size); ++ci)
	{
		const double c = (fractal_params.c.size > 1) ? fractal_params.c.min + (fractal_params.c.max-fractal_params.c.min)/int(fractal_params.c.size-1)*ci : (fractal_params.c.min+fractal_params.c.max)/2;
		fprintf(stderr, "Frame %i of %i. c=%f \n", ci, int(fractal_params.c.size), c);
		render_lyapunov(img[2], fractal_params, *fractal_type, ci);
		
		fprintf(stderr,"\tblurring...\n");
		write_3df_frame(of, img, fractal_params.b.size, fractal_params.a.size);
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
		grid_options = GridOptions::parse(grid_spec);

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
