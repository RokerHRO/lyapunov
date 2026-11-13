#include "rgb.hh"
#include "draw.hxx"

template void draw_hline<RGB>(std::vector<RGB>&, int, int, int, int, int, const RGB&);
template void draw_vline<RGB>(std::vector<RGB>&, int, int, int, int, int, const RGB&);
template void draw_line<RGB>(std::vector<RGB>&, int, int, int, int, int, int, const RGB&);
template void draw_text<RGB>(std::vector<RGB>&, int, int, int, int, const std::string&, const RGB&, const Font& font);
template void draw_shaded_text<RGB>(std::vector<RGB>& img, int width, int height, int x, int y, const std::string& text, const RGB& foreground, const RGB& shadow, const Font& font);
