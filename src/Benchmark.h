/**
 * @file Benchmark.h
 * @brief Headless benchmark mode for the qMandelbrot rendering engine.
 *
 * Provides the command line front end (`--benchmark` and friends) and a
 * runner that renders a fixed view repeatedly without ever showing a window,
 * so the escape-time engine can be timed in isolation from user interaction.
 * The main use is comparing rendering throughput between builds - different
 * compilers, precision modes or iteration limits - on identical input.
 */

#pragma once

#include <QString>
#include <QStringList>
#include "QMandelbrotWidget.h"

/**
 * @struct BenchmarkOptions
 * @brief Configuration of a single benchmark run, populated from the command line.
 *
 * The defaults describe a deep-zoom fixed-point workload: 10 frames of
 * 1280x720 at zoom 2^40 with a 2048 iteration ceiling per pixel. Zoom 2^40 is
 * still inside the range IEEE 754 doubles resolve (the engine switches at
 * 2^44), so the same view can be timed on either precision for a fair
 * comparison between the two engines.
 */
struct BenchmarkOptions {
    bool enabled = false;                                                                  ///< True when --benchmark was given.
    QMandelbrotWidget::Precision precision = QMandelbrotWidget::Precision::FixedPoint128;  ///< Rendering precision under test.
    int64_t pixelIterationLimit = 2048;                                                    ///< Escape-time iteration ceiling per pixel, 0 for auto.
    int32_t imageCount = 10;                                                               ///< Number of frames to render and time.
    int32_t width = 1280;                                                                  ///< Frame width in pixels.
    int32_t height = 720;                                                                  ///< Frame height in pixels.
    int32_t log2Zoom = 40;                                                                 ///< Log2 of the zoom level of the benchmark view.
    bool writeJson = false;                                                                ///< True when --json was given.
    QString jsonPath;                                                                      ///< Destination file for the JSON report.
};

/**
 * @brief Parse the application command line into a BenchmarkOptions.
 *
 * Handles --help and --version internally (both terminate the process, which is
 * what QCommandLineParser::process() does). On a malformed value an error is
 * written to stderr and the function fails without terminating, leaving the
 * caller to decide the exit code.
 *
 * @param arguments Application argument list, as returned by QCoreApplication::arguments().
 * @param options Receives the parsed configuration; untouched fields keep their defaults.
 * @return True when the command line was valid.
 */
bool ParseCommandLine(const QStringList& arguments, BenchmarkOptions& options);

/**
 * @brief Run the benchmark described by @p options and report the results.
 *
 * Renders @c options.imageCount frames offscreen, timing each one, then writes
 * the per-frame timings and their min/max/average/median to the terminal and,
 * when @c options.writeJson is set, to a JSON file.
 *
 * @param options Benchmark configuration.
 * @return Process exit code: 0 on success, non-zero when the JSON report could not be written.
 */
int RunBenchmark(const BenchmarkOptions& options);
