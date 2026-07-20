# lyapunov – a 2D/3D Lyapunov fractal generator
<img src="https://raw.githubusercontent.com/RokerHRO/lyapunov/master/seq-AB.png" width="250" height="250" align="right" />

These programs are in this repo:
## lyapunov
### Requirements
(none, except a libc with getopt() support)
### Compilation

Use `cmake`.

### Running
```
Generate Lyapunov fractal(s) as PPM or 3df to stdout or file(s)

Usage: lyapunov [-W width] [-H height] [-f frames] [-s sequence] [-i max_iter]
     [-t fractal_type ] [-o output_filename_base] [-O format]
   { [-c center_x:center_y:size] | [-x min_x] [-X max_x] [-y min_y] [-Y max_y] } [-z min_z] [-Z max_z]
     [-C value] [-D value] [-E value] [-F value]

        -W width  : width of image in pixel (default: 800)
        -H height : height of image in pixel (default: 800)
        -f frames : number of frames or size in in z direction
        -s seq    : AB or ABC sequence for fractal (default: "BBBABBAAAAAA")
        -i iter   : iterations per pixel (default: 1000)

        -t type    : fractal type (default: plain, run "-t help" to get a list of types)
        -o outfile : basename of the output file (stdout if not given)
                     Can contain template patterns: {s}: sequence, {i}: iter, {f}: frame_nr
        -O format  : either "pnm" for image (series) or "3df" for 3D voxel file (default: "pnm")

        -c cx:cy:x_size : center coordinates and size.
           The x, X, y and Y values will be calculated from them automatically:
        -x min_x  : minimum x value (default: 2.400000)
        -X max_x  : maximum x value (default: 3.999000)
        -y min_y  : minimum y value (default: 3.999000)
        -Y max_y  : maximum y value (default: 2.400000)
        -z min_z  : 3rd dimension: (start) value of z. (selected by "C" in sequence string. default: 2.500000)
        -Z max_z  : maximum value of z (default: 4.000000)
        -C value  : value for 3th dimension (default 2.500000)
        -D value  : value for 4th dimension (default 3.700000)
        -E value  : value for 5th dimension (default 3.000000)
```

<img src="https://raw.githubusercontent.com/RokerHRO/lyapunov/master/seq-ABAAAAB.png" width="250" height="250" align="right" />

#### Fractal types
* `plain`  : just a 2D image or an animation of a linear flight through the 3D cube
* `6sides` : animated cruise through all 6 sides of the ABC cube.
* `circle` : values for C and D form a circle: Center is at `(z_min + z_max)/2`, D.  Radius is `abs(z_max-z_min)/2`. 
             Try this example: `lyapunov -W 320 -H 320 -s ABCDBBCDDA -c 3.48:3.48:1 -z 3.1 -Z 3.6 -D 3.375 -o ani-{s}-{f} -f 320 -t circle`
* `sphere` : a sphere in the ABC cube. D+E can form a circle in animation mode.
             Try this example: `lyapunov -W 1024 -H 512 -t sphere -x 3.2 -X 3.8 -y 3.2 -Y 3.8 -z 3.1 -Z 3.8 -s ABCBBACC -o sphere.ppm`

<img src="https://raw.githubusercontent.com/RokerHRO/lyapunov/master/img/sphere-2222-ABCBBACC-thumb.jpg" width="512" height="256" align="right" />


#### File formats
* ppm : portable pixmap image, or series of images
* 3df : 3D voxel file, suitable for POV-Ray renderer and others

Animations can be viewed via ImageMagick's `animate` command-line tool (if not too big, because all images must fit in RAM!):
```
animate -delay 5 -verbose ani-*.ppm
```
or by converting them to a video file, e.g. with ffmpeg:
```
ffmpeg -i ani-ABCD-%03d.ppm -r 25 ani.mp4
```

## gen-ab-sequences
Prints all non-repetitive AB sequences (that starts with A) to stdout. Sequences starting with B would produce the same images as the sequences with all A's and B's interchanged, except mirrored by the x=y axis.

## gen_zoom
Prints a lot of command lines for the "lyapunov" program to stdout. Running these programs would create a lot of images that shows a lyapunov fractal zoom.

