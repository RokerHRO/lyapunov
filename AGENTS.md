# AGENTS.md

Agent instructions for this repository. Project: Lyapunov fractal generators
(C++20, no dependencies beyond a libc with `getopt()`).

## Programs

- `lyapunov` – the fractal generator (PPM/3df output; `-t help` lists the types)
- `gen-ab-sequences`, `gen_zoom`, `gen_poses` – helper tools (see README)
- `scripts/cube_tour.py` – renders perspective tours of image window poses
  (reads `gen_poses` output; PIL)

## Build

CMake, C++20, any g++-compatible toolchain (on Windows: MSYS2/MinGW, on
Linux: g++/clang++).

- MinGW builds link the dynamic libstdc++ from the toolchain directory:
  that directory must be on PATH **for building AND for running** the
  binaries.
- Keep the build warning-free (unused parameter names are commented out,
  e.g. `/*frame*/`, to document the interface).
- Older MinGW libstdc++ versions lack P0408 operators
  (`string + string_view`, `const char* + string_view`): build such
  strings with explicit `std::string(...)` conversions.

## Code layout and conventions

- `fractal_type.hh/.cc` – all fractal types; `lyapunov.cc` – CLI driver;
  helpers go in the anonymous namespace of `fractal_type.cc`.
- `SimpleFractalType` subclasses override the inline `per_line`/`per_pixel`
  hooks; window-moving types are `FlightPath` subclasses used via `Flight`
  (a static base plus the path).
- Construction goes through the factories `make_fractal_type()` /
  `make_flight_path()` (`std::unique_ptr`, `std::string_view` type, throws
  `std::invalid_argument` for unknown types and invalid parameters – no
  out-parameters).
- Coding style is defined in `.clang-format` (tabs; braces on their own line
  for functions, classes and control statements; left pointer alignment;
  140 columns). The repo is not fully clang-format-clean (the code mixes
  dense operator spacing `a*b` with spaced `a + b`, manually aligned
  help-text/case blocks); keep new code consistent with the surroundings.
- Small focused functions; comments only for non-obvious math (state the
  formula, not the trivial code).

## Adding a new fractal type (checklist)

1. Class in `fractal_type.hh` + implementation in `fractal_type.cc`
   (parameter validation in the constructor/factory, as above)
2. `case` in `make_fractal_type()` (and `make_flight_path()` for
   `FlightPath` types – used by `gen_poses`)
3. Help text in `lyapunov.cc` – it is mirrored 1:1 in the README "Running"
   section; **update both**
4. README: bullet in the "Fractal types" list + an "Examples" section with
   an illustration (GIF: 160x160, 25 fps, ffmpeg `palettegen`/`paletteuse`;
   window tours: `gen_poses` + `scripts/cube_tour.py` -> mp4) and link it
5. New **warp** (window distortion, usable via `-w` with every flight
   type): class implementing `Warp::warp(u, v, ...)` in `fractal_type.hh`,
   working in normalized window coordinates (u = column/b.size horizontal,
   v = line/a.size vertical, both in [0,1]), `case` in `make_warp()`,
   help text + README "Warps" bullet. Warps are applied in
   `Flight::per_pixel` before the pose mapping; they must not change the
   result when their amplitude is 0 (bit-exact).

## Naming of committed media (`img/`)

- fractal animations: `anim-<type>-<sequence>.gif`
- window tours: `<type>-<variant>-tour.mp4` / `<type>-<variant>.mp4`

## Verification workflow

- Render small frames (e.g. 32x32, `-f 4`) and compare bit-exactly
  (hash compare); for math changes cross-check against an independent
  Python implementation before scaling up.
- Before changing `Flight`/`FlightPath`/`Warp` code: rebuild the previous
  binary and keep it (e.g. `lyapunov-baseline.exe` in the temp build
  area); after the change, render the same command lines with both
  binaries and byte-compare (all flight types, single frame and
  animation). Chaotic output amplifies 1-ULP differences, so
  "visually the same" is NOT acceptable.
- `gen_poses` is the source of truth for window poses; the tour rendering
  must stay driven by pose files (no duplicate path math in Python).
- Never write generated output into the repository; use a temp directory.
  Encoded previews (GIF/MP4) need ffmpeg on PATH.
- Windows/PowerShell gotchas: `> file` redirection writes UTF-16 (Python:
  `open(..., encoding="utf-16")`); write multi-line Python to a `.py` file;
  `Select-Object -First` can break native-command pipelines (exit -1, a
  test artifact, not a bug).

## Git

- Commit only when explicitly asked; keep commits small and descriptive.
- The user edits files (especially `README.md`) between turns – re-read a
  file before modifying it, do not assume the last state you saw.

## Known latent bugs (deliberately unfixed)

- `tvals[32][5]`/`tnext[32]` overflow for `OMP_NUM_THREADS > 32`
