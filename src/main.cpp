/**
 * @file main.cpp
 * @brief Application entry point for the qMandelbrot fractal renderer.
 *
 * Parses the command line and then either runs the headless benchmark or
 * initializes the Qt6 application and displays the main window.
 */

#include "pch.h"
#include <QApplication>
#include "Benchmark.h"
#include "QtMainWindow.h"

#ifdef Q_OS_WIN
#include <cstdio>
#include <windows.h>

/**
 * @brief Give this GUI-subsystem process a console to write to.
 *
 * The binary links as a Windows GUI application, so when it is launched from a shell it gets
 * no standard streams at all and anything written to stdout is discarded. Attaching to the
 * console of the launching process routes the benchmark report (and the parser's --help and
 * --version output) to the terminal the command was typed in; when there is no such console
 * one is allocated so the output is visible at all.
 *
 * Redirected streams are left alone: a pipe or a file handed down by the shell
 * (`Mandelbrot.exe --benchmark > run.txt`) is inherited even by a GUI subsystem process and
 * is exactly where the user asked the output to go.
 *
 * A shell does not wait for a GUI application to finish, so its prompt returns immediately
 * and console output arrives interleaved with it.
 */
static void AttachToParentConsole()
{
    // Null std handles are the "GUI app with nowhere to write" case. Anything else - a
    // console, a pipe, a file - is already usable and must not be redirected out from under.
    const HANDLE currentOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (currentOut != nullptr && currentOut != INVALID_HANDLE_VALUE) {
        return;
    }

    if (!AttachConsole(ATTACH_PARENT_PROCESS) && !AllocConsole()) {
        return;
    }

    FILE* stream = nullptr;
    freopen_s(&stream, "CONOUT$", "w", stdout);
    freopen_s(&stream, "CONOUT$", "w", stderr);
}
#endif  // Q_OS_WIN

/**
 * @brief Application entry point.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit code from the benchmark or from the Qt event loop.
 */
int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QCoreApplication::setApplicationName(QStringLiteral("qMandelbrot"));
    QCoreApplication::setApplicationVersion(QStringLiteral("1.0"));

#ifdef Q_OS_WIN
    // Any command line at all means the process was driven from a terminal, and the user is
    // waiting to read what --help, --version or --benchmark have to say.
    if (argc > 1) {
        AttachToParentConsole();
    }
#endif

    BenchmarkOptions options;
    if (!ParseCommandLine(QCoreApplication::arguments(), options)) {
        return 1;
    }

    if (options.enabled) {
        return RunBenchmark(options);
    }

    QtMainWindow w;
    w.show();

    return app.exec();
}
