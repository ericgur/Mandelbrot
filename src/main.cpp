#include "pch.h"
#include <QApplication>
#include "QtMainWindow.h"

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QtMainWindow w;
    w.show();
    return app.exec();
}
