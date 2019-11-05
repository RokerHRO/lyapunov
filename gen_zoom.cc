#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <unistd.h> // for getopt()
#include <vector>

std::string seq = "BBBABBAAAAAA";

double amin = 3.999;
double amax = 2.4;
double bmin = 2.4;
double bmax = 3.999;
int asize = 800;
int bsize = 800;
int nmax = 1000;                /* number of rounds */
int frames = 100;
double zoom = 20;

/* for color generation; somewhat empirical, in order to match Wickerprint's original colors */
double lambda_min = -2.55;
double lambda_max = 0.3959;

int main(int argc, char** argv)
{
	int c;
	while( (c=getopt(argc, argv, "W:H:x:y:X:Y:s:i:F:z:")) != -1 )
	{
		switch(c)
		{
			case 'W': bsize = atoi(optarg); break;
			case 'H': asize = atoi(optarg); break;
			case 'x': bmin = atof(optarg);  break;
			case 'y': amin = atof(optarg);  break;
			case 'X': bmax = atof(optarg);  break;
			case 'Y': amax = atof(optarg);  break;
			case 's': seq = optarg;         break;
			case 'i': nmax = atoi(optarg);  break; 
			case 'F': frames = atoi(optarg); break;
			case 'z': zoom = atof(optarg);  break;
			case '?':
				fprintf(stderr, 
					"Generate commands for a Lyapunov fractal zoom.\n\n"
					"Usage: %s [-W width] [-H height] [-x min_x] [-X max_x] [-y min_y] [-Y max_y] \n"
					          "[-s sequence] [-i max_iter] [-z max_zoomlevel] \n\n"
					"\t-W width  : width of image in pixel (default: %i)\n"
					"\t-H height : height of image in pixel (default: %i)\n"
					"\t-x min_x  : minimum x value (default: %f)\n"
					"\t-X max_x  : maximum x value (default: %f)\n"
					"\t-y min_y  : minimum y value (default: %f)\n"
					"\t-Y max_y  : maximum y value (default: %f)\n"
					"\t-s seq    : AB sequence for fractal (default: \"%s\")\n"
					"\t-i iter   : iterations per pixel (default: %i)\n"
					"\t-F frames : number of frames for zoom flight (default %i)\n"
					"\t-z zoom   : maximum zoom level at end of zoom flight (default: %f)\n"
					, argv[0], bsize, asize,
					bmin, bmax, amin, amax, seq.c_str(), nmax,
					frames, zoom
					);
					return 1;
		}
	}
	
	const double center_b = (bmin + bmax)/2; // horizontal direction
	const double center_a = (amin + amax)/2; // vertical direction
	
	const double db = (bmax - center_b);
	const double da = (amax - center_a);
	
	const double scale = std::pow(1.0/zoom, 1.0/frames);
	for(int f = 0; f<frames; ++f)
	{
		const double scl = std::pow(scale, double(f));
		printf("# Frame #%d:  scale=%f: \n", f, scl );
		printf("./lyapunov -W %d -H %d -s %s -i %i ", bsize, asize, seq.c_str(), nmax );
		printf("-x %10.8f -X %10.8f -y %10.8f -Y %10.8f ",
			(center_b-db*scl), (center_b+db*scl),
			(center_a-da*scl), (center_a+da*scl)
			);
		printf(" | convert -verbose ppm:- zoom-%04i.png\n", f);
	}
}
