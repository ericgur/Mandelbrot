/**
 * @file Benchmark.cpp
 * @brief Implementation of the qMandelbrot headless benchmark mode.
 *
 * Contains the command line definition, the render/timing loop and the two
 * report writers (terminal and JSON).
 */

#include "pch.h"
#include "Benchmark.h"

#include <QCommandLineParser>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTextStream>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

//
// The benchmark always renders the same point of the complex plane, a spot on the boundary in
// the seahorse valley that stays detailed at every zoom level, so the amount of work per frame
// is a function of the configuration alone and two runs are comparable. The coordinates are
// given as decimal strings rather than double literals: past zoom 2^50 or so a double no longer
// has the digits to name the center, and fixed_point128 parses the full precision from text.
//
static constexpr const char* BENCHMARK_CENTER_X = "-0.743643887037158704752191506114774";
static constexpr const char* BENCHMARK_CENTER_Y = "0.131825904205311970493132056385139";

/// @brief File the JSON report goes to when --json is given without a path.
static constexpr const char* DEFAULT_JSON_PATH = "benchmark.json";

/**
 * @struct BenchmarkStatistics
 * @brief Summary of a set of per-frame render times, all in milliseconds.
 */
struct BenchmarkStatistics {
    double min = 0.0;      ///< Fastest frame.
    double max = 0.0;      ///< Slowest frame.
    double average = 0.0;  ///< Arithmetic mean over all frames.
    double median = 0.0;   ///< Middle value, the mean of the two middle values for an even count.
    double total = 0.0;    ///< Sum of all frame times.
};

/**
 * @brief Describe the compiler that built this binary.
 *
 * Reported alongside the timings because telling two toolchains apart is the
 * main reason to run the benchmark. clang-cl is distinguished from a plain
 * clang by the MSVC compatibility macro it also defines.
 *
 * @return Human readable compiler name and version.
 */
[[nodiscard]] static QString CompilerDescription()
{
#if defined(__clang__)
    // __clang_version__ carries the upstream repository URL and commit hash behind the version
    // itself, which is more than a report line needs; the three version macros are just the number.
    const QString clangVersion = QString("%1.%2.%3").arg(__clang_major__).arg(__clang_minor__).arg(__clang_patchlevel__);
#if defined(_MSC_VER)
    return QString("clang-cl ") + clangVersion;
#else
    return QString("clang ") + clangVersion;
#endif
#elif defined(_MSC_VER)
    // _MSC_FULL_VER packs the version as MMmmBBBBB: two digits of major, two of minor and
    // five of build, so 195136256 is 19.51.36256.
    const int32_t version = _MSC_FULL_VER;
    return QString("MSVC %1.%2.%3").arg(version / 10000000).arg((version / 100000) % 100).arg(version % 100000);
#elif defined(__GNUC__)
    return QString("GCC " __VERSION__);
#else
    return QString("unknown");
#endif
}

/**
 * @brief Report the number of threads the render loop will be handed by OpenMP.
 * @return Thread count, or 1 when the build has no OpenMP support.
 */
[[nodiscard]] static int32_t OpenMpThreadCount()
{
#ifdef _OPENMP
    return static_cast<int32_t>(omp_get_max_threads());
#else
    return 1;
#endif
}

/**
 * @brief Convert a precision mode to the token the command line uses for it.
 * @param precision Precision mode.
 * @return Lower case name of the mode.
 */
[[nodiscard]] static QString PrecisionName(QMandelbrotWidget::Precision precision)
{
    switch (precision) {
    case QMandelbrotWidget::Precision::Auto:
        return QString("auto");
    case QMandelbrotWidget::Precision::Double:
        return QString("double");
    case QMandelbrotWidget::Precision::FixedPoint128:
        return QString("fp128");
    case QMandelbrotWidget::Precision::Perturbation:
        return QString("perturbation");
    }

    return QString("unknown");
}

/**
 * @brief Parse a --precision value.
 * @param text Value as given on the command line, case insensitive.
 * @param precision Receives the parsed mode, untouched on failure.
 * @return True when the value named a known mode.
 */
static bool ParsePrecision(const QString& text, QMandelbrotWidget::Precision& precision)
{
    const QString value = text.trimmed().toLower();

    if (value == "auto") {
        precision = QMandelbrotWidget::Precision::Auto;
    } else if (value == "double") {
        precision = QMandelbrotWidget::Precision::Double;
    } else if (value == "fp128" || value == "fixed_point128" || value == "fixed") {
        precision = QMandelbrotWidget::Precision::FixedPoint128;
    } else if (value == "perturbation" || value == "perturb") {
        precision = QMandelbrotWidget::Precision::Perturbation;
    } else {
        return false;
    }

    return true;
}

/**
 * @brief Parse a WIDTHxHEIGHT value, as taken by --size.
 * @param text Value as given on the command line, for example "1920x1080".
 * @param width Receives the width, untouched on failure.
 * @param height Receives the height, untouched on failure.
 * @return True when both dimensions parsed as positive integers.
 */
static bool ParseSize(const QString& text, int32_t& width, int32_t& height)
{
    const QStringList parts = text.trimmed().toLower().split(QChar('x'), Qt::SkipEmptyParts);
    if (parts.size() != 2) {
        return false;
    }

    bool widthOk = false;
    bool heightOk = false;
    const int32_t parsedWidth = parts[0].toInt(&widthOk);
    const int32_t parsedHeight = parts[1].toInt(&heightOk);
    if (!widthOk || !heightOk || parsedWidth <= 0 || parsedHeight <= 0) {
        return false;
    }

    width = parsedWidth;
    height = parsedHeight;

    return true;
}

/**
 * @brief Read an integer option, reporting a range violation on stderr.
 * @param parser Parser holding the parsed command line.
 * @param name Long name of the option.
 * @param minimum Smallest accepted value.
 * @param maximum Largest accepted value.
 * @param result Receives the value, untouched unless the option was given and valid.
 * @return True when the option was absent, or present with a value inside the range.
 */
static bool ReadIntOption(const QCommandLineParser& parser, const QString& name, int64_t minimum, int64_t maximum, int64_t& result)
{
    if (!parser.isSet(name)) {
        return true;
    }

    bool ok = false;
    const int64_t value = parser.value(name).toLongLong(&ok);
    if (!ok || value < minimum || value > maximum) {
        QTextStream(stderr) << "Error: --" << name << " expects an integer in [" << minimum << ", " << maximum << "], got '" << parser.value(name) << "'\n";
        return false;
    }

    result = value;

    return true;
}

/**
 * @brief Give a bare --json argument the default report path.
 *
 * QCommandLineParser has no notion of an option whose value is optional, so a bare --json
 * would be rejected as a missing value. Rewriting it to --json=<default> before parsing keeps
 * both spellings working: --json on its own, and --json <file> to pick the destination.
 *
 * @param arguments Application argument list.
 * @return The list with a bare --json expanded to carry the default path.
 */
[[nodiscard]] static QStringList ExpandBareJsonOption(const QStringList& arguments)
{
    QStringList expanded = arguments;

    for (qsizetype i = 0; i < expanded.size(); ++i) {
        if (expanded[i] != "--json" && expanded[i] != "-json") {
            continue;
        }

        // A path follows only when there is a next argument and it is not an option itself.
        const bool hasValue = (i + 1 < expanded.size()) && !expanded[i + 1].startsWith(QChar('-'));
        if (!hasValue) {
            expanded[i] = QString("--json=") + DEFAULT_JSON_PATH;
        }
    }

    return expanded;
}

bool ParseCommandLine(const QStringList& arguments, BenchmarkOptions& options)
{
    QCommandLineParser parser;
    parser.setApplicationDescription(QString("qMandelbrot - interactive Mandelbrot and Julia set explorer.\n"
                                             "Run without arguments for the GUI, or with --benchmark to time the rendering engine."));
    parser.addHelpOption();
    parser.addVersionOption();

    const QCommandLineOption benchmarkOption(QString("benchmark"), QString("Render offscreen and report timings instead of showing the GUI."));
    const QCommandLineOption precisionOption(QString("precision"), QString("Rendering precision: auto, double, fp128 or perturbation. Default: fp128."),
                                             QString("mode"));
    const QCommandLineOption iterationLimitOption(
        QString("pixel-iteration-limit"), QString("Escape-time iteration ceiling per pixel, or 0 to scale it with zoom. Default: 2048."), QString("count"));
    const QCommandLineOption iterationsOption(QString("iterations"), QString("Number of images to render and time. Default: 10."), QString("count"));
    const QCommandLineOption sizeOption(QString("size"), QString("Image size as WIDTHxHEIGHT. Default: 1280x720."), QString("WxH"));
    const QCommandLineOption zoomOption(QString("zoom"), QString("Log2 of the zoom level of the benchmark view, 0 to 113. Default: 40."), QString("log2"));
    const QCommandLineOption jsonOption(QString("json"), QString("Also write the results to a JSON file. The path is optional and defaults to benchmark.json."),
                                        QString("file"));

    parser.addOption(benchmarkOption);
    parser.addOption(precisionOption);
    parser.addOption(iterationLimitOption);
    parser.addOption(iterationsOption);
    parser.addOption(sizeOption);
    parser.addOption(zoomOption);
    parser.addOption(jsonOption);

    // process() handles --help and --version and exits on a malformed command line.
    parser.process(ExpandBareJsonOption(arguments));

    options.enabled = parser.isSet(benchmarkOption);

    if (parser.isSet(precisionOption) && !ParsePrecision(parser.value(precisionOption), options.precision)) {
        QTextStream(stderr) << "Error: --precision expects auto, double, fp128 or perturbation, got '" << parser.value(precisionOption) << "'\n";
        return false;
    }

    if (parser.isSet(sizeOption) && !ParseSize(parser.value(sizeOption), options.width, options.height)) {
        QTextStream(stderr) << "Error: --size expects WIDTHxHEIGHT, got '" << parser.value(sizeOption) << "'\n";
        return false;
    }

    int64_t iterationLimit = options.pixelIterationLimit;
    int64_t imageCount = options.imageCount;
    int64_t log2Zoom = options.log2Zoom;

    // 0 iterations per pixel is the widget's "scale the limit with the zoom level" mode.
    if (!ReadIntOption(parser, QString("pixel-iteration-limit"), 0, INT32_MAX, iterationLimit) ||
        !ReadIntOption(parser, QString("iterations"), 1, INT32_MAX, imageCount) ||
        !ReadIntOption(parser, QString("zoom"), static_cast<int64_t>(QMandelbrotWidget::logMinZoom), static_cast<int64_t>(QMandelbrotWidget::logMaxZoom),
                       log2Zoom)) {
        return false;
    }

    options.pixelIterationLimit = iterationLimit;
    options.imageCount = static_cast<int32_t>(imageCount);
    options.log2Zoom = static_cast<int32_t>(log2Zoom);

    options.writeJson = parser.isSet(jsonOption);
    options.jsonPath = options.writeJson ? parser.value(jsonOption) : QString();

    return true;
}

/**
 * @brief Summarize a set of per-frame render times.
 * @param timings Per-frame times in milliseconds; may be given in any order.
 * @return The min, max, average, median and total of @p timings, all zero for an empty input.
 */
[[nodiscard]] static BenchmarkStatistics Summarize(const QList<double>& timings)
{
    BenchmarkStatistics stats;
    if (timings.isEmpty()) {
        return stats;
    }

    QList<double> sorted = timings;
    std::sort(sorted.begin(), sorted.end());

    const qsizetype count = sorted.size();
    stats.min = sorted.first();
    stats.max = sorted.last();
    stats.total = std::accumulate(sorted.cbegin(), sorted.cend(), 0.0);
    stats.average = stats.total / static_cast<double>(count);
    stats.median = (count % 2 == 1) ? sorted[count / 2] : (sorted[count / 2 - 1] + sorted[count / 2]) / 2.0;

    return stats;
}

/**
 * @brief Round a millisecond figure to the three decimals the reports carry.
 * @param value Value to round.
 * @return @p value rounded to microsecond resolution.
 */
[[nodiscard]] static double RoundMs(double value)
{
    return std::round(value * 1000.0) / 1000.0;
}

/**
 * @brief Write the JSON report.
 * @param options Benchmark configuration, echoed into the report.
 * @param timings Per-frame render times in milliseconds, in render order.
 * @param stats Summary of @p timings.
 * @param path Destination file.
 * @return True when the file was written.
 */
static bool WriteJsonReport(const BenchmarkOptions& options, const QList<double>& timings, const BenchmarkStatistics& stats, const QString& path)
{
#ifdef NDEBUG
    const QString buildType = QString("release");
#else
    const QString buildType = QString("debug");
#endif

    QJsonObject configuration;
    configuration[QString("compiler")] = CompilerDescription();
    configuration[QString("buildType")] = buildType;
    configuration[QString("precision")] = PrecisionName(options.precision);
    configuration[QString("width")] = options.width;
    configuration[QString("height")] = options.height;
    configuration[QString("log2Zoom")] = options.log2Zoom;
    configuration[QString("pixelIterationLimit")] = static_cast<qint64>(options.pixelIterationLimit);
    configuration[QString("imageCount")] = options.imageCount;
    configuration[QString("openMpThreads")] = OpenMpThreadCount();
    configuration[QString("centerX")] = QString(BENCHMARK_CENTER_X);
    configuration[QString("centerY")] = QString(BENCHMARK_CENTER_Y);

    QJsonArray timingsArray;
    for (double timing : timings) {
        timingsArray.append(RoundMs(timing));
    }

    QJsonObject statistics;
    statistics[QString("min")] = RoundMs(stats.min);
    statistics[QString("max")] = RoundMs(stats.max);
    statistics[QString("average")] = RoundMs(stats.average);
    statistics[QString("median")] = RoundMs(stats.median);
    statistics[QString("total")] = RoundMs(stats.total);

    QJsonObject root;
    root[QString("configuration")] = configuration;
    root[QString("timingsMs")] = timingsArray;
    root[QString("statisticsMs")] = statistics;

    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
        QTextStream(stderr) << "Error: could not write '" << path << "': " << file.errorString() << "\n";
        return false;
    }

    file.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    file.close();

    return true;
}

int RunBenchmark(const BenchmarkOptions& options)
{
    QTextStream out(stdout);

    // The widget is never shown: resize() alone gives it the geometry the renderer works from,
    // and renderOffscreen() drives the full pipeline without a paint event.
    QMandelbrotWidget widget;
    widget.resize(options.width, options.height);
    widget.setMaximumIterations(options.pixelIterationLimit);
    widget.setPrecision(options.precision);
    widget.setView(fp128_t(BENCHMARK_CENTER_X), fp128_t(BENCHMARK_CENTER_Y), options.log2Zoom);

    out << "qMandelbrot benchmark\n";
    out << "  compiler              : " << CompilerDescription() << "\n";
    out << "  precision             : " << PrecisionName(options.precision) << "\n";
    out << "  image size            : " << options.width << " x " << options.height << "\n";
    out << "  zoom                  : 2^" << options.log2Zoom << "\n";
    out << "  pixel iteration limit : " << (options.pixelIterationLimit == 0 ? QString("auto") : QString::number(options.pixelIterationLimit)) << "\n";
    out << "  images                : " << options.imageCount << "\n";
    out << "  OpenMP threads        : " << OpenMpThreadCount() << "\n";
    out << "\n";
    out.flush();

    QList<double> timings;
    timings.reserve(options.imageCount);

    for (int32_t i = 0; i < options.imageCount; ++i) {
        // Every frame has to recompute the fractal from scratch, otherwise the cached
        // iteration buffer would make all but the first frame a palette conversion.
        widget.invalidate();

        const double elapsedMs = widget.renderOffscreen();
        timings.append(elapsedMs);

        out << QString("  image %1 / %2 : %3 ms\n").arg(i + 1, 4).arg(options.imageCount).arg(elapsedMs, 10, 'f', 3);
        out.flush();
    }

    const BenchmarkStatistics stats = Summarize(timings);

    out << "\n";
    out << "Results (ms)\n";
    out << "  min     : " << QString::number(stats.min, 'f', 3) << "\n";
    out << "  max     : " << QString::number(stats.max, 'f', 3) << "\n";
    out << "  average : " << QString::number(stats.average, 'f', 3) << "\n";
    out << "  median  : " << QString::number(stats.median, 'f', 3) << "\n";
    out << "  total   : " << QString::number(stats.total, 'f', 3) << "\n";
    out.flush();

    if (options.writeJson) {
        if (!WriteJsonReport(options, timings, stats, options.jsonPath)) {
            return 1;
        }

        out << "\nJSON report written to " << options.jsonPath << "\n";
        out.flush();
    }

    return 0;
}
