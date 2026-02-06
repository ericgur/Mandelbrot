# qMandelbrot

A high-performance Mandelbrot and Julia set fractal renderer built with Qt6 and modern C++20. Features dual-precision rendering with a custom 128-bit fixed-point arithmetic library, enabling extreme zoom levels up to 2^113.

## Features

- **Dual rendering precision** -- standard 64-bit double and custom 128-bit fixed-point arithmetic with automatic switching based on zoom depth
- **Mandelbrot and Julia sets** -- switch between fractal types; Julia set includes 10 curated presets and custom constant input
- **Multiple color palettes** -- Grey, Gradient, Vivid (HSV rainbow), and Histogram-equalized coloring with smooth iteration interpolation
- **Palette animation** -- real-time color cycling with configurable timer
- **OpenMP parallelization** -- multi-threaded scanline rendering with dynamic scheduling, toggleable at runtime
- **Image export** -- save rendered fractals as PNG at 1920x1080, 2560x1440, or 3840x2160
- **Interactive navigation** -- mouse click to zoom/pan, keyboard arrow keys for panning, +/- for zoom
- **Cross-platform** -- builds on Windows (MSVC, Clang/LLVM), macOS, and Linux

## Requirements

- **C++20** compiler (MSVC 2019+, Clang 12+, or GCC 9+)
- **CMake** 3.25 or later
- **Qt6** (Core, Gui, Widgets)
- **OpenMP** runtime
- **Ninja** build system (used by the provided scripts)

## Building

### Configure and build

**Linux / macOS:**
```bash
./configure.sh
./build.sh
```

**Windows:**
```bat
configure.bat
build.bat
```

The executable is placed in the `bin/` directory as `qMandelbrot` (or `qMandelbrot.exe` on Windows).

### Other build scripts

| Script | Description |
|--------|-------------|
| `clean.sh` / `clean.bat` | Delete the `build/` directory |
| `rebuild.sh` / `rebuild.bat` | Clean then build |
| `install.sh` / `install.bat` | Run CMake install target |

### Environment setup

- **macOS:** Source `qt_env_macos.sh` to configure Qt paths before building.
- **Windows (LLVM):** Run `qt_env_llvm_windows.bat` to set up the Clang/LLVM Qt environment.

## Usage

### Mouse controls

| Action | Effect |
|--------|--------|
| Left click | Zoom in 2x centered on cursor |
| Ctrl + Left click | Zoom in 4x |
| Ctrl + Shift + Left click | Zoom in 8x |
| Right click | Zoom out 2x |
| Ctrl + Right click | Zoom out 4x |
| Ctrl + Shift + Right click | Zoom out 8x |
| Middle click | Reset to default view |

### Keyboard controls

| Key | Effect |
|-----|--------|
| Arrow keys | Pan view (5% of viewport per press) |
| `+` | Zoom in 2x (centered) |
| `-` | Zoom out 2x (centered) |

### Menu options

- **File** -- Save rendered image at preset resolutions (PNG), Exit
- **View**
  - **Precision** -- Auto (switches at zoom > 2^44), Double, FixedPoint128
  - **Set Type** -- Mandelbrot or Julia
  - **Iterations** -- Auto (scales with zoom level), or fixed: 128, 192, 256, 384, 512, 768, 1024, 1536, 2048
  - **Reset Zoom** -- return to initial view
  - **Animate Palette** -- toggle color cycling animation
  - **Julia Set Options** -- open dialog to select presets or enter custom constants
  - **OpenMP** -- toggle parallel rendering

## Project Structure

```
.
├── CMakeLists.txt              # Root CMake configuration
├── src/
│   ├── CMakeLists.txt          # Source-level build config
│   ├── main.cpp                # Application entry point
│   ├── pch.h / pch.cpp         # Precompiled header
│   ├── QtMainWindow.h/.cpp     # Main window (menus, actions, status bar)
│   ├── QtMainWindow.ui         # Qt Designer UI for main window
│   ├── QMandelbrotWidget.h/.cpp# Core fractal rendering engine
│   ├── QJuliaSetOptions.h/.cpp # Julia set constant configuration dialog
│   ├── QJuliaSetOptions.ui     # Qt Designer UI for Julia options
│   ├── fixed_point128.h        # 128-bit fixed-point arithmetic library
│   ├── fixed_point128_shared.h # Shared definitions for fixed-point types
│   ├── fixed_point128.natvis   # Visual Studio debugger visualizer for fp128
│   ├── Mandelbrot.rc           # Windows resource file (icon)
│   ├── resource.h              # Windows resource IDs
│   └── res/icons/              # Application icons (.ico, .icns, .png)
├── bin/                        # Compiled binaries
├── build/                      # CMake build artifacts (generated)
├── *.sh / *.bat                # Build, configure, clean, install scripts
├── Mandelbrot.sln              # Visual Studio solution
└── LICENSE                     # MIT License
```

## Source Files

### `main.cpp`

Application entry point. Creates a `QApplication` and shows the `QtMainWindow`.

### `QtMainWindow` (`QtMainWindow.h` / `QtMainWindow.cpp`)

The main application window derived from `QMainWindow`. Manages the menu bar, action groups, keyboard input, and status bar display. Connects UI actions to the rendering widget for precision selection, iteration count, set type switching, and image saving. Displays render statistics (time, zoom level, resolution, iterations) in the status bar after each frame.

### `QMandelbrotWidget` (`QMandelbrotWidget.h` / `QMandelbrotWidget.cpp`)

The core rendering engine derived from `QWidget`. Owns the complete fractal computation and display pipeline:

- **Fractal computation** -- Implements the escape-time algorithm for both Mandelbrot (`Z(0) = 0, C = pixel`) and Julia (`Z(0) = pixel, C = constant`) sets. Iterates `Z = Z^2 + C` until `|Z| > 2` or the iteration limit is reached.
- **Dual precision** -- `DrawImageDouble()` uses IEEE 754 doubles; `DrawImageFixedPoint128()` uses the custom `fp128_t` type (8 integer bits, 120 fractional bits) for deep zoom. In Auto mode, switches to fixed-point when zoom exceeds 2^44.
- **Color mapping** -- Supports four palette types (Grey, Gradient, Vivid, Histogram). Uses continuous smoothing via `mu = iter - log(log(|Z|)) / log(2)` and alpha-blended interpolation between adjacent colors.
- **Auto-iterations** -- Scales the iteration limit linearly with `log2(zoom)` from 128 (at 1x) to 2500 (at 2^113 zoom).
- **OpenMP** -- Parallelizes the main rendering loop per scanline with `#pragma omp parallel for schedule(dynamic)`.

### `QJuliaSetOptions` (`QJuliaSetOptions.h` / `QJuliaSetOptions.cpp`)

A dialog for configuring Julia set parameters. Provides 10 preset complex constants (e.g., `0.285 + 0.01i`, `-0.8 + 0.156i`) via a combo box, plus manual real/imaginary input validated to the range [-2, 2]. Supports auto-apply mode for real-time constant changes.

### `fixed_point128.h` / `fixed_point128_shared.h`

A header-only 128-bit fixed-point arithmetic library. The template class `fixed_point128<I>` parameterizes the number of integer bits (the project uses `I=8`, giving 120 fractional bits). Provides:

- Arithmetic operators: `+`, `-`, `*`, `/`, shifts, comparisons
- Math functions: `sqrt`, `sin`, `cos`, `tan`, `atan`, `atan2`, `exp`, `log`, `log2`, `pow`, and more
- Utility functions: `fabs`, `floor`, `ceil`, `trunc`, `round`, `fmod`, `copysign`, `hypot`
- Platform-specific intrinsics for MSVC (`_mulx_u64`, `_addcarryx_u64`, `_udiv128`) and GCC (`__int128`, `__builtin_*`)
- Visual Studio debugger visualization via `fixed_point128.natvis`

**Acknowledgements** (from the header):
- `div_32bit` derived from *Hacker's Delight* 2nd Edition by Henry S. Warren Jr.
- Logarithm functions from Dan Moulding's [log2fix](https://github.com/dmoulding/log2fix)
- Square root based on *Math Toolkit for Real Time Programming* by Jack W. Crenshaw

### `pch.h` / `pch.cpp`

Precompiled header including `<memory>` and `<QApplication>` to speed up builds.

## Rendering Algorithm

The Mandelbrot set is computed by iterating the complex recurrence:

```
Z(0) = 0,  C = (x, y)          -- Mandelbrot
Z(0) = (x, y),  C = constant   -- Julia

Z(i) = Z(i-1)^2 + C
```

Iteration continues until `|Z| > 2` (escape radius) or the maximum iteration count is reached. Points that never escape are colored black (set interior).

**Smooth coloring** uses the normalized iteration count:

```
mu = iter + 1 - log(log(|Z|)) / log(2)
```

This fractional iteration value indexes into the color table with linear interpolation between adjacent entries, eliminating visible banding.

## License

MIT License. Copyright (c) 2022 Eric Gur. See [LICENSE](LICENSE) for details.
