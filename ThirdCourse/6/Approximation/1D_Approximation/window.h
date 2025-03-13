#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  const char *f_name;
  double a;
  double b;
  int n;
  double (*f) (double);

  // Перечисление для типов аппроксимации
  enum ApproximationType {
    SPLINE,
    CHEBYSHEV,
    BOTH
  };

  ApproximationType currentApproximation; // Текущий тип аппроксимации

public:
  Window (QWidget *parent);

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);

public slots:
  void change_func ();
  void toggle_approximation(); // Метод для переключения типа аппроксимации

protected:
  void paintEvent (QPaintEvent *event);
  void keyPressEvent(QKeyEvent *event) override; // Обработчик нажатий клавиш
};

#endif