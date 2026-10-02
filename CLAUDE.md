# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**qMandelbrot** is a Qt6 GUI application for interactive Mandelbrot/Julia set exploration with dual-precision rendering (IEEE 754 double and custom 128-bit fixed-point arithmetic for extreme zoom depths).

## Build Commands

**Prerequisites:** CMake 3.25+, Qt6, OpenMP, Ninja

```bash
# First-time setup
./configure.sh          # Linux/macOS
configure.bat           # Windows

# Build (default: RelWithDebInfo)
./build.sh              # Linux/macOS
./build.sh Debug        # specify config
build.bat               # Windows

# Other
./rebuild.sh            # clean + build
./install.sh            # install to /bin
./clean.sh              # remove build/ directory

# Qt environment (must be sourced before configure/build)
source qt_env_macos.sh           # macOS
qt_env_llvm_windows.bat          # Windows (LLVM/Clang toolchain)
```

Output binary: `bin/qMandelbrot` (or `bin/qMandelbrot.exe` on Windows)

## Code Formatting

```bash
./apply_format.sh       # Run clang-format on all src/ files
```

Configuration is in `.clang-format`. No automated test suite — testing is done manually through the GUI.

## Architecture

### Rendering Pipeline

The app renders fractals via an escape-time algorithm:
- **Mandelbrot:** `Z(0)=0, C=(x,y)`, iterate `Z = Z² + C` until `|Z| > 2`
- **Julia:** `Z(0)=(x,y), C=constant`, same iteration

Smooth coloring: `mu = iter + 1 - log(log(|Z|)) / log(2)` — eliminates banding via fractional iteration interpolation between palette entries.

### Dual-Precision Strategy (`QMandelbrotWidget`)

- `CalcIterationsDouble()` — IEEE 754 `double`, up to zoom ~2⁴⁴
- `CalcIterationsFP128()` — custom `fixed_point128<fp128IntBits>`, currently 3 integer bits + sign + 124 fraction bits, up to zoom 2^`logMaxZoom` (`F - 10` = 2¹¹⁴, where pixels of a 3840 wide image are still 1 LSB apart)
- **Fewer integer bits = deeper zoom.** `fp128IntBits = 3` holds [-8, 8): view center clamped to ±`maxCenterMagnitude` (2), pixel offsets to ±5, iterates stay below 6.83. The escaping modulus (up to 46.6) wraps, so the escape test `Bounded()` also checks |u|, |v| < 2, as one OR of high QWORDs (measured free; four signed compares cost 7.5%)
- Pixel coordinates come from `PixelCoordinate()`: center + halfWidth * (2i - count) / width, with the ratio divided exactly in fp128, so each is good to a few LSBs. Never step from the edge by a rounded `dx` (its error times the pixel index wrecks images past 2^100), and don't form the ratio in `double` (off by 2^-53 of the half width, which flips the counts of chaotic boundary pixels)
- Auto mode switches precision when zoom exceeds 2⁴⁴

### Key Classes

| Class | File | Role |
|---|---|---|
| `QMandelbrotWidget` | `src/QMandelbrotWidget.h/.cpp` | Core rendering engine; escape-time algo, coloring, mouse/keyboard input, image export |
| `QtMainWindow` | `src/QtMainWindow.h/.cpp` | Main window; menu bar, status bar stats, routes actions to widget |
| `QJuliaSetOptions` | `src/QJuliaSetOptions.h/.cpp` | Dialog for Julia set constant selection (10 presets + manual input) |
| `Complex128` | `src/QMandelbrotWidget.h` | Complex number with `fp128_t` parts; the Julia constant everywhere (`std::complex` is only specified for built-in floats) |
| `Favorite` | `src/Favorites.h/.cpp` | Saved location (center, log2 zoom, iteration limit, set type, Julia constant); QSettings load/save |
| `QFavoritesDialog` | `src/QFavoritesDialog.h/.cpp` | Modal favorites editor (list + detail form; edits a copy, committed on OK) |
| (functions) | `src/Fp128Text.h/.cpp` | `fp128_t` <-> decimal text: validation patterns, parsing, shortest round-trip formatting |
| `fixed_point128<I>` | `src/fixed_point128.h` | Header-only 128-bit fixed-point arithmetic with full math function suite |

### `QMandelbrotWidget` Internals

- **OpenMP parallelization:** `#pragma omp parallel for schedule(dynamic)` per scanline
- **Color palettes:** Grey, Gradient, Vivid, Histogram-equalized
- **View:** `_centerX`, `_centerY` and the zoom level are the whole view state. Zoom is always a power of 2, held as one `int32_t _logZoomLevel` clamped to [`logMinZoom`, `logMaxZoom`]; `ViewHalfWidth()` is `2.5 >> level`, exact in fp128, and the vertical extent follows from the image aspect ratio. A zoom the clamp turns into no change leaves the view untouched
- **Auto-iterations:** Scales iteration limit with `log2(zoom)` from 128 (1x) to 2500 (2^`logMaxZoom`)
- **Color animation:** `QChronoTimer` drives palette cycling
- **Lazy redraw:** `_fractalDataValid` and `_colorTableValid` flags gate recomputation of iteration data and color LUT independently. `invalidate(false)` leaves the color table flag alone; only the color table builders may mark it valid, or a stale table is read past its end after a limit change
- **Mouse:** Left click = zoom 2x in; Right click = zoom 2x out; Middle = reset. Ctrl multiplies by 2x, Ctrl+Shift by 4x
- **Keyboard:** Arrow keys pan 5%; +/- zoom 2x; F8 saves the view as a favorite (handled by `QtMainWindow`)
- **Export:** PNG at 1920×1080, 2560×1440, or 3840×2160

### Favorites

- Stored with `QSettings` (native format; organization `ericgur` set in `main.cpp`, required for QSettings to read or write at all)
- Coordinates and Julia constants are stored as decimal text, never `double`, so they keep fp128 precision; `ToDecimalText()` writes the shortest text that parses back to the same `fp128_t`
- Iteration limit is stored as Auto (`auto_iterations`, written as `auto`) or a fixed value; restoring goes through `QtMainWindow::SelectIterationLimit()` so the Iterations menu and slider stay in sync
- First run (no `favorites/size` key) returns `DefaultFavorites()`, ten famous locations, each with a fixed iteration limit picked by measuring unresolved pixels; an emptied list is stored with size 0 and stays empty
- Applying a favorite must switch set type and Julia constant **before** `setView()`, because `setSetType()` and `setJuliaConstant()` both reset the view

### `fixed_point128<I>` Library

Header-only library in `src/fixed_point128.h`. Supports all standard arithmetic operators and math functions (`sqrt`, `sin`, `cos`, `tan`, `atan2`, `exp`, `log`, `log2`, `pow`, etc.). Uses platform intrinsics:
- MSVC: `_mulx_u64`, `_addcarryx_u64`, `_udiv128`
- GCC/Clang: `__int128`, `__builtin_*`

Debugger visualization: `src/fixed_point128.natvis` (Visual Studio)

## Build Configuration

- **Standard:** C++20 (C17 for C files)
- **Modules:** Qt6 Core, Gui, Widgets
- **Parallelism:** OpenMP
- **Optimizations:** `-flto`, `-march=native`, warnings as errors for format security
- **Windows (MSVC):** AVX2 instruction set

## Code Style Rules

- Document using Doxygen-style comments; class headers require verbose documentation.
- **PascalCase:** global functions, private/protected class methods.
- **camelCase:** public methods, struct methods/members, local variables, and any function starting with `q` (Qt style).
- **Private/protected data members:** underscore prefix + camelCase (e.g., `_dataMember`).
- No newline before opening curly brace — single space only (applies to `if`, `for`, `while`, `try`, etc.).
- All control blocks require curly braces, including single-line bodies with macros or function calls.
- Final `return` statement must be on its own line.
- All non-error/non-bool return values must be marked `[[nodiscard]]`.
- 4-space indentation; no tabs.
