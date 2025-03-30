#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>
#include <cmath>
#include "chebyshov_approximation.hpp"
#include "spline_approximation.hpp"

#define PI M_PI

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id = 0;
  const char *f_name = nullptr;
  double a = -1;
  double b = 1;
  int n = 0;
  int p = 0;
  bool first = true;
  double (*f)(double) = nullptr;

  enum ApproximationType {
    SPLINE,
    CHEBYSHEV,
    BOTH,
    ERRORS
  };

  ApproximationType currentApproximation = SPLINE;

  struct ApproximationData {
    double *alpha = nullptr;
    double *c = nullptr;
    double *A = nullptr;
    double *g = nullptr;
    double *g2 = nullptr;
    double *z = nullptr;
    double* x = nullptr;
    double *F = nullptr;
    double *dF = nullptr;
    bool chebyshevCalculated = false;
    bool splineCalculated = false;

    ~ApproximationData() { clear(); }

    void clear() {
      delete[] alpha; alpha = nullptr;
      delete[] c; c = nullptr;
      delete[] A; A = nullptr;
      delete[] g; g = nullptr;
      delete[] g2; g2 = nullptr;
      delete[] z; z = nullptr;
      delete[] x; x = nullptr;
      delete[] F; F = nullptr;
      delete[] dF; dF = nullptr;
      chebyshevCalculated = false;
      splineCalculated = false;
    }

    void allocate(int n) {
      clear();
      alpha = new double[n];
      c = new double[4*n];
      A = new double[3*n];
      g = new double[n];
      g2 = new double[n];
      z = new double[n];
      x = new double [n];
      F = new double[n];
      dF = new double[n];
    }
  };

  ApproximationData approxData;

  void calculateChebyshevApproximation(double max_y, double min_y);
  void calculateSplineApproximation(double max_y, double min_y);
  void clearApproximationData();

public:
  Window(QWidget *parent);
  QSize minimumSizeHint() const;
  QSize sizeHint() const;
  int parse_command_line(int argc, char *argv[]);
  QPointF l2g(double x_loc, double y_loc, double y_min, double y_max);

public slots:
  void change_func();
  void update_function();
  void toggle_approximation();
  void zoom_in();
  void zoom_out();
  void increase_points();
  void decrease_points();
  void point_down();
  void point_up();

protected:
  void paintEvent(QPaintEvent *event);
  void keyPressEvent(QKeyEvent *event) override;
};

#endif