#pragma once

#include <QMainWindow>
#include "ui_QtMainWindow.h"
#include "QJuliaSetOptions.h"

class QAction;
class QAbstractButton;
class QMandelbrotWidget;

class QtMainWindow : public QMainWindow
{
    Q_OBJECT
public:
    explicit QtMainWindow(QWidget* parent = nullptr);
    ~QtMainWindow();

private slots:
    void onActionSaveImage();
    void onActionPrecision();
    void onActionIterations();
    void onActionSetType();
    void onActionJuliaOptions();
    void onRenderDone(FrameStats stats);

private:
    void createActions();

    Ui_QtMainWindow ui;
    QMandelbrotWidget* m_centralWidget;
    QVector<QAction*> iterActions;
    QActionGroup* setTypeGroup = nullptr;
    QActionGroup* iterGroup = nullptr;
    QJuliaSetOptions* juliaOptionsDialog = nullptr;
};
