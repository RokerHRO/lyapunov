# lyapunov – a 2d/3D Lyapunov fractal generator

Thes programs are in this repo:
## lyapunov
### Requirements
(none, except a libc with getopt() support)
### Compilation
c++ -Wall -O3 -std=c++14 -o lyapunov lyapunov.cc
### Running
```
Usage: ./lyapunov [-W width] [-H height] [-z value] [-s sequence] [-i max_iter]
   { [-c center_x:center_y:size] | [-x min_x] [-X max_x] [-y min_y] [-Y max_y] }

	-W width  : width of image in pixel (default: 800)
	-H height : height of image in pixel (default: 800)
	-F frames : number of frames in z direction (if specified a 3df file is generaed, not a PPM)
	-z value  : 3rd dimension: (start) value of z. (selected by "C" in sequence string. default: 2.500000)
	-Z max_z  : maximum value of z (used for 3df files) 
	-s seq    : AB or ABC sequence for fractal (default: "BBBABBAAAAAA")
	-i iter   : iterations per pixel (default: 1000)

	-c cx:cy:x_size : center coordinates and size.
	   The following values will be calculated from them automatically:
	-x min_x  : minimum x value (default: 2.400000)
	-X max_x  : maximum x value (default: 3.999000)
	-y min_y  : minimum y value (default: 3.999000)
	-Y max_y  : maximum y value (default: 2.400000)
```
The files (PPM image for 2D fractals or Pov-Ray's 3df file format for 3D fractals) are written to stdout.

## gen-ab-sequences
Prints all non-repetitive AB sequences (that starts with A) to stdout. Sequences starting with B would produce the same images as the sequences with all A's and B's interchanged, except mirrored by the x=y axis.

## gen_zoom
Prints a lot of command lines for the "lyapunov" program to stdout. Running these programs would create a lot of images that shows a lyapunov fractal zoom.
