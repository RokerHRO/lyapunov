#include "render.hxx"

template void render_lyapunov<RGB>(std::vector<RGB>&, const FractalParams&, const FractalType&, unsigned frame);
template void render_lyapunov<uint8_t>(std::vector<uint8_t>&, const FractalParams&, const FractalType&, unsigned frame);
