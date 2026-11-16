#ifndef RENDER_HH
#define RENDER_HH

#include <cstddef>
#include <vector>

#include "fractal_type.hh"

// public interface: render one frame of the Lyapunov fractal into a row-major
// pixel buffer of size_t(width) * height pixels (ai = row, bi = column)
template<class PXL>
void render_lyapunov(std::vector<PXL>& img, const FractalParams& params, const FractalType& ftype, unsigned frame);

#endif // RENDER_HH
