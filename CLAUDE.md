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
- `CalcIterationsFP128()` — custom `fixed_point128<8>` (8 int bits, 120 fractional bits), up to zoom ~2¹¹³
- Auto mode switches precision when zoom exceeds 2⁴⁴

### Key Classes

| Class | File | Role |
|---|---|---|
| `QMandelbrotWidget` | `src/QMandelbrotWidget.h/.cpp` | Core rendering engine; escape-time algo, coloring, mouse/keyboard input, image export |
| `QtMainWindow` | `src/QtMainWindow.h/.cpp` | Main window; menu bar, status bar stats, routes actions to widget |
| `QJuliaSetOptions` | `src/QJuliaSetOptions.h/.cpp` | Dialog for Julia set constant selection (10 presets + manual input) |
| `fixed_point128<I>` | `src/fixed_point128.h` | Header-only 128-bit fixed-point arithmetic with full math function suite |

### `QMandelbrotWidget` Internals

- **OpenMP parallelization:** `#pragma omp parallel for schedule(dynamic)` per scanline
- **Color palettes:** Grey, Gradient, Vivid, Histogram-equalized
- **Auto-iterations:** Scales iteration limit with `log2(zoom)` from 128 (1x) to 2500 (2¹¹³)
- **Color animation:** `QChronoTimer` drives palette cycling
- **Lazy redraw:** `_fractalDataValid` and `_colorTableValid` flags gate recomputation of iteration data and color LUT independently
- **Mouse:** Left click = zoom 2x in; Right click = zoom 2x out; Middle = reset. Ctrl multiplies by 2x, Ctrl+Shift by 4x
- **Keyboard:** Arrow keys pan 5%; +/- zoom 2x
- **Export:** PNG at 1920×1080, 2560×1440, or 3840×2160

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
