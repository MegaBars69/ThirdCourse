#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id = 0;
  const char *f_name= nullptr;
  double a = -1;
  double b = 1;
  int n = 0; 
  int p = 0;
  bool first = true;
  double (*f) (double) = nullptr;

  // Перечисление для типов аппроксимации
  enum ApproximationType {
    SPLINE,
    CHEBYSHEV,
    BOTH,
    ERRORS
  };

  ApproximationType currentApproximation = SPLINE; // Текущий тип аппроксимации

public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);

public slots:
  void change_func ();
  void update_function();
  void toggle_approximation(); // Метод для переключения типа аппроксимации
  void zoom_in();  // Уменьшение отрезка [a, b] в 2 раза
  void zoom_out(); // Увеличение отрезка [a, b] в 2 раза
  void increase_points(); // Увеличение количества точек в 2 раза
  void decrease_points(); // Уменьшение количества точек в 2 раза
  void point_down();
  void point_up();
protected:
  void paintEvent (QPaintEvent *event);
  void keyPressEvent(QKeyEvent *event) override; // Обработчик нажатий клавиш
};

#endif