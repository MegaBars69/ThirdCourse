// surfacewindow.h
#ifndef SURFACEWINDOW_H
#define SURFACEWINDOW_H
#define EPSILON 1e-15
#include "algorithm.hpp"
#include "thread_function.hpp"
#include "initialize_matrix.hpp"

#include <QWidget>
#include <QVector>
#include <QVector3D>
#include <QMatrix4x4>
#include <QPainter>
#include <QWheelEvent>
#include <QMouseEvent>

enum CalculationStatus
{
    UNDEF,
    CALCULATING,
    CALCULATED
};

class SurfaceWindow : public QWidget {
    Q_OBJECT
public:
    explicit SurfaceWindow(QWidget *parent = nullptr);
    void calculateSurface();

    QSize minimumSizeHint() const override { return QSize(100, 100); }
    QSize sizeHint() const override { return QSize(1000, 1000); }
    int parse_command_line(int argc, char *argv[]);
    void CalcNorm();
public slots:
    void zoomIn() { zoom *= 1.1; update(); }
    void zoomOut() { zoom *= 0.9; update(); }
    void rotateLeft() { rotationY -= 10; update(); }
    void rotateRight() { rotationY += 10; update(); }
    void change_func();
    void update_function();
    void toggle_approximation();
    void zoom_in();
    void zoom_out();
    void increase_points();
    void decrease_points();
    void increase_draw_points();
    void decrease_draw_points();
    void point_down();
    void point_up();
    void reset_view();  
    void drawCoordinateSystem(QPainter &painter);
    void drawAxis(QPainter &painter, 
                 const QVector3D &start, 
                 const QVector3D &end,
                 const QString &label);
    void drawCube(QPainter &painter);
protected:
    void paintEvent(QPaintEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void clearApproximationData();
    void keyPressEvent(QKeyEvent *event) override;
    void ApproximateFunction();
    void closeEvent(QCloseEvent *event) override;
private:
    double a = -1, b = 1, c = -1, d = 1, eps = 1e-10;
    int nx = 5, ny = 5, func_id = 0, max_it = 100, p=1, its = 0; 
    double r1 = -1, r2 = -1, r3 = -1, r4 = -1, t1 = 0.0, t2 = 0.0;
    int mx = 5, my = 5, point = 0;
    double (*f)(double,double) = nullptr;
    const char *f_name = nullptr;
    bool* threads_working = nullptr;

    bool first = true;
    double m_minZ = -1;
    double m_maxZ = 1;
    double norm = 1;
    pthread_mutex_t p_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t p_cond = PTHREAD_COND_INITIALIZER;
    struct Triangle 
    {
        QVector3D points[3];
        double depth;
    };

    enum CurrentGraph 
    {
        FUNC,
        APPROX,
        ERRORS
    };

    CurrentGraph currentFunc = FUNC;
    QVector<QVector3D> vertices;
    QVector<Triangle> triangles;
    QMatrix4x4 projection;
    double zoom = 1.33;
    double rotationX = -68.5;
    double rotationY = -0.5;
    QPoint lastMousePos;

    QVector3D project(const QVector3D &point) const;
    void updateProjection();

    struct ApproximationData 
    {
        Args* aA = nullptr;
        double* A = nullptr, *B = nullptr, *x = nullptr, *u = nullptr, *v = nullptr, *r = nullptr;
        int *I = nullptr;
        int len_msr = 1;
        int N = 1;
        

        CalculationStatus calc_status = UNDEF;
        
        ~ApproximationData() { clear(); }
    
        void clear();
        void allocate(int nx, int ny, int p);
        
      };
      ApproximationData approxData;
      
};

#endif // SURFACEWINDOW_H