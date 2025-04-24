// main.cpp
#include <QApplication>
#include <QMainWindow>
#include <QMenuBar>
#include <QVBoxLayout>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include "surfacewindow.h"
#include <fenv.h>

int main(int argc, char *argv[]) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    QApplication app(argc, argv);

    QMainWindow *mainWindow = new QMainWindow;
    SurfaceWindow *surface = new SurfaceWindow(mainWindow);
    if (surface->parse_command_line(argc, argv)) {
        QMessageBox::warning(0, "Wrong input arguments!", "Wrong input arguments!");
        return -1;
    }

    // Настройка меню
    QMenuBar *menuBar = new QMenuBar(mainWindow);
    QMenu *viewMenu = menuBar->addMenu("View");

    QAction *zoomIn = viewMenu->addAction("Zoom In");
    QAction *zoomOut = viewMenu->addAction("Zoom Out");
    QAction *rotateLeft = viewMenu->addAction("Rotate Left");
    QAction *rotateRight = viewMenu->addAction("Rotate Right");
    QAction *action;

    // Connect number key actions (0-7)
    action = menuBar->addAction("&Change function (0)", surface, SLOT(change_func()));
    action->setShortcut(QString("0"));

    action = menuBar->addAction("&Toggle approximation (1)", surface, SLOT(toggle_approximation()));
    action->setShortcut(QString("1"));

    action = menuBar->addAction("&Zoom in (2)", surface, SLOT(zoom_in()));
    action->setShortcut(QString("2"));

    action = menuBar->addAction("&Zoom out (3)", surface, SLOT(zoom_out()));
    action->setShortcut(QString("3"));

    action = menuBar->addAction("&Increase points (4)", surface, SLOT(increase_points()));
    action->setShortcut(QString("4"));

    action = menuBar->addAction("&Decrease points (5)", surface, SLOT(decrease_points()));
    action->setShortcut(QString("5"));

    action = menuBar->addAction("&Point up (6)", surface, SLOT(point_up()));
    action->setShortcut(QString("6"));

    action = menuBar->addAction("&Point down (7)", surface, SLOT(point_down()));
    action->setShortcut(QString("7"));

    action = menuBar->addAction("&(mx,my)++ (8)", surface, SLOT(increase_draw_points()));
    action->setShortcut(QString("8"));

    action = menuBar->addAction("&(mx,my)-- (9)", surface, SLOT(decrease_draw_points()));
    action->setShortcut(QString("9"));
    
    action = menuBar->addAction("&Reset View (r)", surface, SLOT(reset_view()));
    action->setShortcut(QString("r"));
   
    action = menuBar->addAction("&Exit", mainWindow, SLOT(close()));
    action->setShortcut(QString("Ctrl+X"));

    QObject::connect(zoomIn, &QAction::triggered, surface, &SurfaceWindow::zoomIn);
    QObject::connect(zoomOut, &QAction::triggered, surface, &SurfaceWindow::zoomOut);
    QObject::connect(rotateLeft, &QAction::triggered, surface, &SurfaceWindow::rotateLeft);
    QObject::connect(rotateRight, &QAction::triggered, surface, &SurfaceWindow::rotateRight);

    menuBar->setMaximumHeight(30);
    mainWindow->setMenuBar(menuBar);
    mainWindow->setCentralWidget(surface);
    surface->setWindowTitle("Graph");
    
    // Инициализация поверхности
    //surface->calculateSurface();

    mainWindow->show();
    app.exec();
    //delete surface;
    return 0;
}