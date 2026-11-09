#include "rgb.hh"
#include "draw.hxx"

template void draw_hline<RGB>(std::vector<RGB>&, const int, const int, const int, const int, const int, const RGB&);
template void draw_vline<RGB>(std::vector<RGB>&, const int, const int, const int, const int, const int, const RGB&);
template void draw_line<RGB>(std::vector<RGB>&, const int, const int, int, int, const int, const int, const RGB&);
template void draw_text<RGB>(std::vector<RGB>&, const int, const int, const int, const int, const std::string&, const RGB&, const RGB&, const int);
