/**
 * @file main.cpp
 * @brief Application entry point for the qMandelbrot fractal renderer.
 *
 * Initializes the Qt6 application and displays the main window.
 */

#include "pch.h"
#include <QApplication>
#include "QtMainWindow.h"

/**
 * @brief Application entry point.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return Exit code from the Qt event loop.
 */
int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QtMainWindow w;
    w.show();
    return app.exec();
}
