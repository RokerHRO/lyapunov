#include "../draw.hh"

#include "ter-16.c"
#include "ter-20.c"
#include "ter-24.c"
#include "ter-32.c"

const std::map<std::string, const Font*> fonts = {
	{"16", &ter_16},
	{"20", &ter_20},
	{"24", &ter_24},
	{"32", &ter_32},
};
