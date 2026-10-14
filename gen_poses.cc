// Dumps the image window poses of a "Flight" fractal type (6sides, wander)
// to stdout or a file. The poses are in the unit ABC cube, exactly the ones
// used by the fractal render, so the visualization script scripts/cube_tour.py
// does not need its own copy of the path math.
//
// Usage: gen_poses -t <6sides|wander> [-p params] [-f frames] [-o outfile]
//
// Output, one line per item (doubles with 17 significant digits, i.e. the
// full round-trip precision):
//   POSE <frame> corner(3) u_axis(3) v_axis(3)
// plus the type-specific auxiliary data of FlightPath::dumpPathData():
//   CURVE <x> <y> <z>      (wander: the closed path of the window centers)
//   KEYFRAME <x> <y> <z>   (wander: the random keyframe centers)
//   HINGE <frame> a(3) b(3) (6sides, ease 0/1: the hinge edge per frame)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h> // for getopt()
#include <vector>

#include "fractal_type.hh"

static std::vector<double> parse_params(const char* p)
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
	int a;
	std::string type;
	std::vector<double> params;
	int frames = 256;
	std::string outfile;

	while( (a=getopt(argc, argv, "t:p:f:o:")) != -1 )
	{
		switch(a)
		{
			case 't': type = optarg; break;
			case 'p': params = parse_params(optarg); break;
			case 'f': frames = atoi(optarg); break;
			case 'o': outfile = optarg; break;
			default:
				fprintf(stderr,
					"Dumps the window poses of a flight fractal type for scripts/cube_tour.py.\n\n"
					"Usage: gen_poses -t <6sides|wander> [-p params] [-f frames] [-o outfile]\n");
				return 2;
		}
	}

	if (frames < 1)
	{
		fprintf(stderr, "The number of frames must be >= 1.\n");
		return 2;
	}

	std::unique_ptr<FlightPath> path;
	try
	{
		path = make_flight_path(type, params);
	}
	catch (const std::invalid_argument& e)
	{
		fprintf(stderr, "%s\n", e.what());
		return 2;
	}

	std::ofstream of;
	if (!outfile.empty())
	{
		of.open(outfile);
		if (!of)
		{
			fprintf(stderr, "Cannot open output file \"%s\"\n", outfile.c_str());
			return 2;
		}
	}
	std::ostream& os = (!outfile.empty()) ? of : std::cout;
	os << std::setprecision(17);

	os << "# gen_poses type=" << type << " frames=" << frames << "\n";
	for (int frame = 0; frame < frames; ++frame)
	{
		const WindowPose p = path->poseAt(unsigned(frame), unsigned(frames));
		os << "POSE " << frame << ' '
			<< p.corner.x << ' ' << p.corner.y << ' ' << p.corner.z << ' '
			<< p.u_axis.x << ' ' << p.u_axis.y << ' ' << p.u_axis.z << ' '
			<< p.v_axis.x << ' ' << p.v_axis.y << ' ' << p.v_axis.z << '\n';
	}
	path->dumpPathData(os, unsigned(frames));
	return 0;
}
