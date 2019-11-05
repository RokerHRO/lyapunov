#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <unistd.h> // for getopt()
#include <vector>

#define CLAMP(x) (((x) > 255) ? (unsigned char)255u : ((x) < 0) ? (unsigned char)0u : (unsigned char)(x))

//int seq_length = 12;
//int seq[12] = {1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0};

std::string center;
std::string seq = "BBBABBAAAAAA";

double cmin = 2.5; // We go to 3rd dimension! \o/
double cmax = 4.0; // for 3df files

double amin = 3.999;
double amax = 2.4;
double bmin = 2.4;
double bmax = 3.999;

int asize = 800;
int bsize = 800;
int csize = 0;  // value >0 lets create 3df file instead of PPM

int nmax = 1000;                /* number of rounds */

/* for color generation; somewhat empirical, in order to match Wickerprint's original colors */
double lambda_min = -2.55;
double lambda_max = 0.3959;

struct RGB
{
	unsigned char r,g,b;
};


template<class PXL>
PXL lambda2pxl(double lambda);

template<>
RGB lambda2pxl<RGB>(double lambda)
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
uint8_t lambda2pxl<uint8_t>(double lambda)
{
	return CLAMP( 127 - lambda*100 );
}


// extension: 3rd dimension: value c is fixed, selected by 'C' in sequence string
template<class PXL>
void lyapunov(std::vector<std::vector<PXL>>& img, int max_iter, const std::string& sequence, double c)
{

  double min_x = +HUGE_VAL;
  double max_x = -HUGE_VAL;
  double min_lambda = +HUGE_VAL;
  double max_lambda = -HUGE_VAL;

  const int height = img.size();
  const int width  = img[0].size();

#pragma omp parallel for schedule(dynamic)
  for (int ai = 0; ai <height; ++ai)
  {
    const double a = amin + (amax-amin)/height*(ai+0.5);
    double sum_log_deriv, prod_deriv, x;
    fprintf(stderr, "\r%4d of %4d   ", ai, height);
    for (int bi = 0; bi < width; ++bi)
    {
      const double b = bmin + (bmax-bmin)/width*(bi+0.5);
      x = 0.5;
      /* one round without derivating, so that the value 0.5 is avoided */
      for (unsigned m = 0; m < sequence.length(); m++)
      {
        const char ch = sequence[m];
        const double r = ch=='B' ? b : (ch=='A' ? a : c);
        x = r*x*(1-x);
      }
      sum_log_deriv = 0;
      for (int n = 0; n < max_iter; n++)
      {
        prod_deriv = 1;
        for (unsigned m = 0; m < sequence.length(); m++)
        {
          const char ch = sequence[m];
          const double r = ch=='B' ? b : (ch=='A' ? a : c);
          /* avoid computing too many logarithms. One every round is acceptable. */
          prod_deriv *= r*(1-2*x); 
          x = r*x*(1-x);
        }
        sum_log_deriv += log(fabs(prod_deriv));
      }
      const double lambda = sum_log_deriv / (nmax*sequence.length());
      
      if(x>max_x) max_x = x;
      if(x<min_x) min_x = x;
      if(lambda>max_lambda) max_lambda = lambda;
      if(lambda<min_lambda) min_lambda = lambda;
      
      img[ai][bi] = lambda2pxl<PXL>(lambda);
    }
  }
  
  fprintf(stderr, " x=%f ... %f ; lambda=%f ... %f \n",
      min_x, max_x, min_lambda, max_lambda );
  
}


void make_ppm()
{
	std::vector<std::vector<RGB>> img;
	img.resize(asize);
	for(auto& row : img)
	{
		row.resize(bsize);
	}
	
	// PPM header
	printf("P6\n");
	fflush(stdout);
	lyapunov(img, nmax, seq, cmin);
	printf("# Lyapunov: max_iter=%d  seq='%s'  a=%f ... %f  b=%f ... %f  c=%f  "
		"%d %d\n255\n",
		nmax, seq.c_str(), amin, amax, bmin, bmax, cmin,
		bsize, asize
		);
	
	for(const auto& row : img)
	{
		fwrite(row.data(), sizeof(RGB), row.size(), stdout);
	}
	fflush(stdout);
}

inline
uint8_t get(const std::vector<std::vector<std::vector<uint8_t>>>& img, int x, int y, int z)
{
//	if( unsigned(z)>=img.size() ) return 0;
	if( unsigned(y)>=img[z].size() ) return 0;
	if( unsigned(x)>=img[z][y].size() ) return 0;
	return img[z][y][x];
}


void make_3df()
{
	// 3df header
	printf("%c%c%c%c%c%c",
		bsize >> 8, bsize & 255,
		asize >> 8, asize & 255,
		csize >> 8, csize & 255
	);
	fflush(stdout);
	
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
		lyapunov(img[2], nmax, seq, c);
		
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
			fputc(u8, stdout);
		}
		fflush(stdout);
		img[0].swap(img[1]);
		img[1].swap(img[2]);
	}
}


int main(int argc, char** argv)
{
	int a;
	while( (a=getopt(argc, argv, "W:H:F:z:Z:c:x:y:X:Y:s:i:")) != -1 )
	{
		switch(a)
		{
			case 'W': bsize = atoi(optarg); break;
			case 'H': asize = atoi(optarg); break;
			case 'F': csize = atoi(optarg); break;
			case 'z': cmin  = atof(optarg); break;
			case 'Z': cmax  = atof(optarg); break;
			case 'c': center = optarg;      break;
			case 'x': bmin = atof(optarg);  break;
			case 'y': amin = atof(optarg);  break;
			case 'X': bmax = atof(optarg);  break;
			case 'Y': amax = atof(optarg);  break;
			case 's': seq = optarg;         break;
			case 'i': nmax = atoi(optarg);  break; 
			case '?':
				fprintf(stderr, 
					"Generate a Lyapunov fractal as PPM or 3df to stdout.\n\n"
					"Usage: %s [-W width] [-H height] [-z value] [-s sequence] [-i max_iter]\n"
					"   { [-c center_x:center_y:size] | [-x min_x] [-X max_x] [-y min_y] [-Y max_y] }\n"
					"\n"
					"\t-W width  : width of image in pixel (default: %i)\n"
					"\t-H height : height of image in pixel (default: %i)\n"
					"\t-F frames : number of frames in z direction (if specified a 3df file is generaed, not a PPM)\n"
					"\t-z value  : 3rd dimension: (start) value of z. (selected by \"C\" in sequence string. default: %f)\n"
					"\t-Z max_z  : maximum value of z (used for 3df files) \n"
					"\t-s seq    : AB or ABC sequence for fractal (default: \"%s\")\n"
					"\t-i iter   : iterations per pixel (default: %i)\n"
					"\n"
					"\t-c cx:cy:x_size : center coordinates and size.\n"
					"\t   The following values will be calculated from them automatically:\n"
					"\t-x min_x  : minimum x value (default: %f)\n"
					"\t-X max_x  : maximum x value (default: %f)\n"
					"\t-y min_y  : minimum y value (default: %f)\n"
					"\t-Y max_y  : maximum y value (default: %f)\n"
					, argv[0], bsize, asize, cmin, seq.c_str(), nmax,
					bmin, bmax, amin, amax
					);
					return 1;
		}
	}
	
	if(center.size())
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
		
		fprintf(stderr, "Interval calculated from \"-c %s\": min_x=%11.9f  max_x=%11.9f  min_y=%11.9f  max_y=%11.9f\n",
			center.c_str(), bmin, bmax, amin, amax);
	}
	
	if(csize)
	{
		make_3df();
	}else{
		make_ppm();
	}
	
	fflush(stdout);
	fprintf(stderr, "\n");
}
