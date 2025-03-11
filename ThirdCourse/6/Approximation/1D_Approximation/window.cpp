
#include <QPainter>
#include <stdio.h>
#include <iostream>
#include <cmath>
#define PI 3.14159265358980

#include "approximation.hpp"
#include "window.h"

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))

static
double f_0 (double x)
{
  x= x*x;
  return 1;
}

static
double f_1 (double x)
{
  return x;
}

static
double f_2 (double x)
{
  return x*x;
}

static
double f_3 (double x)
{
  return x*x*x;
}

static
double f_4 (double x)
{
  return x*x*x*x;
}

static
double f_5 (double x)
{
  return pow(2.718281828459045,x);
}
static
double f_6 (double x)
{
  return (1/(25*x*x+1));
}
Window::Window (QWidget *parent)
  : QWidget (parent)
{
  a = DEFAULT_A;
  b = DEFAULT_B;
  n = DEFAULT_N;

  func_id = 0;

  change_func ();
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

int Window::parse_command_line (int argc, char *argv[])
{
  if (argc == 1)
    return 0;

  if (argc == 2)
    return -1;

  if (   sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || (argc > 3 && sscanf (argv[3], "%d", &n) != 1)
      || n <= 0)
    return -2;

  return 0;
}

/// change current function for drawing
void Window::change_func ()
{
  func_id = (func_id + 1) % 7;

  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;
        break;
      case 2:
        f_name = "f (x) = x*x";
        f = f_2;
        break;
      case 3:
        f_name = "f (x) = x * x * x";
        f = f_3;
        break;
      case 4:
        f_name = "f (x) = x*x*x*x";
        f = f_4;
        break;
      case 5:
        f_name = "f (x) = e^x";
        f = f_5;
        break;
      case 6:
        f_name = "f (x) = 1/(25*x*x+1)";
        f = f_6;
        break;
    }
  update ();
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max)
{
  double x_gl = (x_loc - a) / (b - a) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{  
  QPainter painter (this);
  int M = 1080;
  double x1, x2, y1, y2;
  double max_y, min_y;
  //double step = 0;
  double xm = 0;
  double delta_y, delta_x = (b - a) / M;
  QPen pen_black(Qt::black, 0, Qt::SolidLine); 
  QPen pen_red(Qt::red, 0, Qt::SolidLine); 
  QPen pen_blue(Qt::blue, 3, Qt::SolidLine); 
  
  painter.setPen (pen_blue); 
  double *alpha = new double[n];
  double *g = new double[n];
  double *g2 = new double[n];
  double *z = new double[n];
  double*F = new double[n];
  memset(alpha,0,n*sizeof(double));
  memset(g,0,n*sizeof(double));
  memset(g2,0,n*sizeof(double));
  memset(z,0,n*sizeof(double));

  //step =(b-a)/(n-1); 
  for (int m = 0; m < n; m++)
  {
    xm = (a+b)/2 + (b-a)*cos(PI*(2*m+1)/(2*n));
    F[m] = f(xm);
  }
  ChebyshovAproximation(n, F, alpha, g, g2, z);

  // calculate min and max for current function
  max_y = min_y = 0;
  for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
    {
      y1 = ChebyshovValue(x1, a,b,n,alpha);
      //y1 = f(x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  // draw approximated line for graph
  x1 = a;
  y1 = ChebyshovValue(x1, a,b,n,alpha);  
  //y1 = f(x1);
  for (x2 = x1 + delta_x; x2 - b < 1.e-6; x2 += delta_x) 
    {
      y2 = ChebyshovValue(x2, a,b,n,alpha);
      //y2 = f(x2);
      // local coords are converted to draw coords
      painter.drawLine (L2G(x1, y1), L2G(x2, y2));

      x1 = x2, y1 = y2;
    }
  x2 = b;
  //y2 = f(x2);
  y2 = ChebyshovValue(x2, a,b,n,alpha);
  painter.drawLine (L2G(x1, y1), L2G(x2, y2));
 
  // draw axis
  painter.setPen (pen_black);
  painter.drawLine (L2G(a, 0), L2G(b, 0));
  painter.drawLine (L2G(0, min_y), L2G(0, max_y));

  // render function name
  painter.setPen ("black");
  painter.drawText (0, 20, f_name);

  delete[] alpha;
  delete[] z;
  delete[] g; 
  delete[] g2;
  delete[] F;

}
