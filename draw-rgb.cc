#include "rgb.hh"
#include "draw.hxx"

template void draw_hline<RGB>(std::vector<RGB>&, int, int, int, int, int, const RGB&);
template void draw_vline<RGB>(std::vector<RGB>&, int, int, int, int, int, const RGB&);
template void draw_line<RGB>(std::vector<RGB>&, int, int, int, int, int, int, const RGB&);
template void draw_text<RGB>(std::vector<RGB>&, int, int, int, int, const std::string&, const RGB&, const Font& font);
