#include "fractal_params.hh"
#include "string_case.hh"

std::string FractalParams::frame_nr_to_string(const unsigned frame_nr) const
{
	if (c.size == 0)
		return {};
	char buf[16] = {};
	const unsigned digits = snprintf(nullptr, 0, "%u", c.size-1);
	snprintf(buf, sizeof(buf), "%0*u", digits, frame_nr);
	return buf;
}

std::string FractalParams::expand_patterns(const std::string tmpl, const unsigned frame_nr) const
{
	std::string t = tmpl;
	t = replace_all(t, "{s}", sequence);
	t = replace_all(t, "{i}", std::to_string(max_iter));
	t = replace_all(t, "{f}", frame_nr_to_string(frame_nr));
	return t;
}
