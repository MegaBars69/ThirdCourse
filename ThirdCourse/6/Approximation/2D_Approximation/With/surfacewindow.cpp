#include "surfacewindow.h"
#include <cmath>
#include <algorithm>
#include <QWheelEvent>
#include <QCloseEvent>
#include <QMouseEvent>
#include <QPainterPath>
#include <QMessageBox>
#include <QPixmap>
#include "algorithm.hpp"
#include "thread_function.hpp"
#include "initialize_matrix.hpp"

static
double f_0(double /*x*/, double)
{
    return 1;
}

static
double f_1(double x, double)
{
    return x;
}

static
double f_2(double /*x*/, double y)
{
    return y;
}

static
double f_3(double x, double y)
{
    return x + y;
}

static
double f_4(double x, double y)
{
    return sqrt(x*x+y*y);
}

static
double f_5(double x, double y)
{
    return x*x+y*y;
}

static
double f_6(double x, double y)
{
    return exp(x*x-y*y);
}
static
double f_7(double x, double y)
{
    return 1/(25*(x*x + y*y) + 1);
}

void SurfaceWindow::change_func()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }
    
    if (!first)
    {
        func_id = (func_id + 1) % 8;
    }
    else
    {
        first = false;
    }
    norm = 0;
    strochka_printed = false;

    currentFunc = FUNC;
    clearApproximationData();
    update_function();
    update();
}

void SurfaceWindow::update_function()
{
    switch (func_id)
    {
        case 0:
            f_name = "f(x, y) = 1";
            f = f_0;
            break;
        case 1:
            f_name = "f(x, y) = x";
            f = f_1;
            break;
        case 2:
            f_name = "f(x, y) = y";
            f = f_2;
            break;
        case 3:
            f_name = "f(x, y) = x + y";
            f = f_3;
            break;
        case 4:
            f_name = "f(x, y) = sqrt(x*x+y*y)";
            f = f_4;
            break;
        case 5:
            f_name = "f(x, y) = x*x+y*y";
            f = f_5;
            break;
        case 6:
            f_name = "f(x, y) = e^(x*x-y*y)";
            f = f_6;
            break;
        case 7:
            f_name = "f(x, y) = 1/(25(x*x + y*y)+1)";
            f = f_7;
            break;
    }
}

void SurfaceWindow::toggle_approximation()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }

    if (currentFunc == FUNC)
        currentFunc = APPROX;
    else if (currentFunc == APPROX)
        currentFunc = ERRORS;
    else
    {
        currentFunc = FUNC;
    }
    strochka_printed = false;
    update();
}


void SurfaceWindow::zoom_in()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }

    double center = (a + b) / 2.0;
    double half_length = (b - a) / 4.0;
    a = center - half_length;
    b = center + half_length;
    point = 0;
    
    center = (c + d) / 2.0;
    half_length = (d - c) / 4.0;
    c = center - half_length;
    d = center + half_length;
    strochka_printed = false;

    clearApproximationData();
    update();
}

void SurfaceWindow::zoom_out()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }

    double center = (a + b) / 2.0;
    double half_length = (b - a);
    a = center - half_length;
    b = center + half_length;

    center = (c + d) / 2.0;
    half_length = (d - c);
    c = center - half_length;
    d = center + half_length;
    strochka_printed = false;

    clearApproximationData();
    update();
}

void SurfaceWindow::increase_points()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }

    nx *= 2;
    ny *= 2;
    strochka_printed = false;

    clearApproximationData();
    update();
}

void SurfaceWindow::decrease_points()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }

    if (nx > 3 && ny > 3)
    {
        nx /= 2;
        ny /= 2;
        strochka_printed = false;
        clearApproximationData();
        update();
    }
}
void SurfaceWindow::increase_draw_points()
{
    mx *= 2;
    my *= 2;
    update();
}

void SurfaceWindow::decrease_draw_points()
{
    if (mx > 3 && my > 3)
    {
        mx /= 2;
        my /= 2;
        update();
    }
}
void SurfaceWindow::point_down()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }
    point--;
    strochka_printed = false;
    clearApproximationData();
    update();
}

void SurfaceWindow::point_up()
{
    if (threads_working[0] == true) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }
    point++;
    strochka_printed = false;

    clearApproximationData();
    update();
}


void SurfaceWindow::keyPressEvent(QKeyEvent *event) {
    // Проверяем статус вычислений
    // Если вычислений нет - обрабатываем клавиши как обычно
    switch(event->key()) {
        case Qt::Key_0: change_func(); break;            
        case Qt::Key_1: toggle_approximation(); break;
        case Qt::Key_2: zoom_in(); break;
        case Qt::Key_3: zoom_out(); break;
        case Qt::Key_4: increase_points(); break;
        case Qt::Key_5: decrease_points(); break;
        case Qt::Key_6: point_up(); break;
        case Qt::Key_7: point_down(); break;
        case Qt::Key_8: increase_draw_points(); break;
        case Qt::Key_9: decrease_draw_points(); break;
        case Qt::Key_R: reset_view(); break;
        case Qt::Key_Left: rotateLeft(); break;   // Вращение влево
        case Qt::Key_Right: rotateRight(); break; // Вращение вправо
        case Qt::Key_Up: rotateClockwise(); break; // Вращение по часовой стрелке
        case Qt::Key_Down: rotateCounterClockwise(); break; // Вращение против часовой стрелки
        default: QWidget::keyPressEvent(event);
    }

}



void SurfaceWindow::rotateClockwise() {
    rotationZ += 10; // Увеличиваем угол вращения по оси Z
    update(); // Обновляем окно
}

void SurfaceWindow::rotateCounterClockwise() {
    rotationZ -= 10; // Уменьшаем угол вращения по оси Z
    update(); // Обновляем окно
}
void SurfaceWindow::reset_view() {
    zoom = 1.33;
    rotationX = -68.5;
    rotationY = -0.5;
    update();
}

void drawSmiley(QPainter &painter, int x, int y, int size)
{
    // Рисуем лицо
    painter.drawEllipse(x, y, size, size);
    /*
    // Рисуем глаза
    painter.setBrush(Qt::black);
    painter.drawEllipse(x + size / 4, y + size / 4, size / 8, size / 8); // Левый глаз
    painter.drawEllipse(x + 3 * size / 4 - size / 8, y + size / 4, size / 8, size / 8); // Правый глаз
       
    // Рисуем улыбку
    QPainterPath smile;
    smile.moveTo(x + size / 4, y + size / 2);
    smile.quadTo(x + size / 2, y + 3 * size / 4, x + 3 * size / 4, y + size / 2);
    painter.drawPath(smile);
    */
}

void SurfaceWindow::CalcNorm()
{
    double xi, yj, z = 0;

    double hx_loc = (b - a) / mx;
    double hy_loc = (d - c) / my;
   
    // Инициализация min/max
    m_minZ = std::numeric_limits<double>::max();
    m_maxZ = std::numeric_limits<double>::lowest();
    
    // Generate vertices
    for (int i = 0; i <= mx; ++i) 
    {
        for (int j = 0; j <= my; ++j) 
        {
            xi = a + i * hx_loc;
            yj = c + j * hy_loc;
            z = f(xi,yj);
            
            // Обновляем min/max
            if (z < m_minZ) m_minZ = z;
            if (z > m_maxZ) m_maxZ = z;
        }
    }
    norm = (fabs(m_maxZ) > fabs(m_minZ) ? fabs(m_maxZ) : fabs(m_minZ));

}

void SurfaceWindow::calculateSurface() 
{
    pthread_mutex_lock(&p_mutex); 
    if(approxData.calc_status == CALCULATING)
    {
        approxData.calc_status = (*threads_working ? CALCULATING : CALCULATED);
    }
    pthread_mutex_unlock(&p_mutex);
    vertices.clear();
    triangles.clear();

    double xi, yj, z = 0;

    double hx_loc = (b - a) / mx;
    double hy_loc = (d - c) / my;
    double hx = (b - a) / nx;
    double hy = (d - c) / ny;

    if(fabs(norm) < EPSILON)
    {
        CalcNorm();
    }
    if(currentFunc != FUNC && approxData.calc_status != CALCULATED)
    {
        approxData.calc_status = CALCULATING;
        ApproximateFunction();
    }
    m_minZ = std::numeric_limits<double>::max();
    m_maxZ = std::numeric_limits<double>::lowest();
    // Generate vertices
    for (int i = 0; i <= mx; ++i)
    {
        for (int j = 0; j <= my; ++j)
        {
            xi = a + i * hx_loc;
            yj = c + j * hy_loc;
            if(currentFunc == FUNC)
            {
                z = f(xi,yj);
            }
            else if(currentFunc == APPROX && approxData.calc_status == CALCULATED && j < my && i < mx)
            {
                z = Pf(approxData.x, xi, yj, a, c, hx, hy, nx, ny);
            }
            else if(currentFunc == ERRORS && approxData.calc_status == CALCULATED && j < my && i < mx)
            {
                z = fabs(f(xi,yj) - Pf(approxData.x, xi, yj, a, c, hx, hy, nx, ny));
            }
            else if(currentFunc == ERRORS && approxData.calc_status == CALCULATED && (j == my || i == mx))
            {
                z = fabs(f(xi,yj) - f(xi, yj));
            }
            else if(currentFunc == APPROX && approxData.calc_status == CALCULATING)
            {
                z = f(xi,yj);
            }
            else if(currentFunc == ERRORS && approxData.calc_status == UNDEF)
            {
                z = f(xi,yj);
            }
            else
            {
                z = f(xi,yj);
            }
            
            if(i < mx && j < my)
            {
                // Обновляем min/max
                if (z < m_minZ) m_minZ = z;
                if (z > m_maxZ) m_maxZ = z;
            }

            vertices.append(QVector3D(xi, yj, z));
        }
    }
    max = fabs(m_maxZ) > fabs(m_minZ) ? fabs(m_maxZ) : fabs(m_minZ);

    // Защита от деления на ноль
    if (qFuzzyCompare(m_minZ, m_maxZ)) {
        m_maxZ = m_minZ + 0.001;
    }

    // Create triangles
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            int idx = i * (my + 1) + j;
            Triangle t1, t2;
            
            t1.points[0] = vertices[idx];
            t1.points[1] = vertices[idx + 1];
            t1.points[2] = vertices[idx + my + 1];
            
            t2.points[0] = vertices[idx + 1];
            t2.points[1] = vertices[idx + my + 2];
            t2.points[2] = vertices[idx + my + 1];
            
            triangles.append(t1);
            triangles.append(t2);
        }
    }
    update();
}
void SurfaceWindow::paintEvent(QPaintEvent *) 
{
    QPainter painter(this);
    
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(), Qt::white);
    calculateSurface();

    update_function();
    updateProjection();
    QPen pen_black(Qt::black, 2, Qt::SolidLine);
    double normalizedZ, delta_y;

    const char* prefix = "max{|Fmax||,|Fmin|} = ";
    char* strochka = new char[strlen(prefix) + 20];
    for (Triangle &tri : triangles) {
        QVector3D center = (tri.points[0] + tri.points[1] + tri.points[2]) / 3;
        tri.depth = (projection * center).z();
    }
    
    std::sort(triangles.begin(), triangles.end(), 
        [](const Triangle &a, const Triangle &b) { return a.depth > b.depth; });
    
    drawCoordinateSystem(painter);

    // Draw triangles
    for (const Triangle &tri : triangles) {
        QPolygonF polygon;
        for (const QVector3D &p : tri.points) {
            QVector3D proj = project(p);
            polygon << QPointF(proj.x(), proj.y());
        }
        
        double avgZ = (tri.points[0].z() + tri.points[1].z() + tri.points[2].z()) / 3;
        normalizedZ = (avgZ - m_minZ) / (m_maxZ - m_minZ);
        normalizedZ = qBound(0.0, normalizedZ, 1.0); 
        // Цвет от синего (min) к красному (max)
        QColor color = QColor::fromHsvF(0.7 * (1.0 - normalizedZ), 0.8, 1.0);
        
        painter.setBrush(color);
        painter.setPen(Qt::NoPen);
        painter.drawPolygon(polygon);
    }
    QFont font2("Arial", 12, QFont::Bold);
    painter.setFont(font2);
    painter.setPen(pen_black);
    painter.drawText(0, 20, f_name);
    
    delta_y = 0.01 * (m_maxZ - m_minZ);
    prefix = "max{|Fmax||,|Fmin|} = ";
    m_minZ -= delta_y;
    m_maxZ += delta_y;
    
    sprintf(strochka, "%s%e", prefix, max);
    if(!strochka_printed && approxData.calc_status != CALCULATING)
    {
        printf("%s\n", strochka);  
        strochka_printed = true;
    }
    painter.drawText(0, 130, strochka);

    prefix = "(nx, ny) = ";
    sprintf(strochka, "(nx, ny) = (%d, %d)", nx, ny);
    painter.drawText(0, 95, strochka);

    // For (mx, my)
    sprintf(strochka, "(mx, my) = (%d, %d)", mx, my);
    painter.drawText(0, 115, strochka);

    prefix = "p = ";
    sprintf(strochka, "%s%d", prefix, point);
    painter.drawText(0, 155, strochka);
    
    painter.setPen(pen_black);

    if (currentFunc == FUNC)
    {
        painter.drawText(0, 45, "FUNCTION");
    }
    else if (currentFunc == APPROX)
    {
        painter.drawText(0, 45, "APPROXIMATION");
    }
    else if (currentFunc == ERRORS)
    {
        painter.setPen(pen_black);
        painter.drawText(0, 45, "ERRORS");
    }
    QPixmap controlsPixmap("mouse.png"); // Путь к изображению
    QSize newSize(100, 100); // Укажите желаемый размер
    QPixmap scaledPixmap = controlsPixmap.scaled(newSize, Qt::KeepAspectRatio, Qt::SmoothTransformation);

    // Отображение в правом верхнем углу
    int xPos = width() - scaledPixmap.width() - 10; 
    int yPos = 10; 
    painter.drawPixmap(xPos, yPos, scaledPixmap); 
    painter.setPen(Qt::black); 
    painter.setFont(QFont("Arial", 12, QFont::Bold)); 
    painter.drawText(xPos-25, yPos + scaledPixmap.height()+ 20, "Mouse Available"); 

    QPixmap arrowsPixmap("arrows.png"); 
    QSize new2Size(125, 125); 
    QPixmap scaledarrowsPixmap = arrowsPixmap.scaled(new2Size, Qt::KeepAspectRatio, Qt::SmoothTransformation);

    xPos = width() - scaledarrowsPixmap.width()+10; 
    yPos = 150; 
    painter.drawPixmap(xPos, yPos, scaledarrowsPixmap);

    painter.setPen(Qt::black); 
    painter.setFont(QFont("Arial", 12, QFont::Bold));
    painter.drawText(xPos-25, yPos + scaledarrowsPixmap.height()+ 20, "Arrows Available");
    delete[] strochka;
}

QVector3D SurfaceWindow::project(const QVector3D &point) const {
    QMatrix4x4 view;
    view.translate(0, 0, -5 * zoom);
    view.rotate(rotationX, 1, 0, 0);
    view.rotate(rotationY, 0, 1, 0);
    view.rotate(rotationZ, 0, 0, 1); 

    QVector3D proj = projection * view * point;
    proj.setX((proj.x() + 1) * width() / 2);
    proj.setY((1 - proj.y()) * height() / 2);
    return proj;
}

void SurfaceWindow::updateProjection()
{
    projection.setToIdentity();
    
    QSize viewport = size();
    double aspect = static_cast<double>(viewport.width()) / viewport.height();
    
    QVector3D sceneSize(b - a, d - c, m_maxZ - m_minZ);
    double maxSceneSize = qMax(sceneSize.x(), qMax(sceneSize.y(), sceneSize.z()));
    
    maxSceneSize *= zoom;
    
    projection.ortho(
        -maxSceneSize * aspect/2, maxSceneSize * aspect/2, 
        -maxSceneSize/2,         maxSceneSize/2,         
        -maxSceneSize*2,          maxSceneSize*2         
    );
    
    projection.translate(
        -((a + b)/2 + c + d)/2, 
        -((c + d)/2 + m_minZ + m_maxZ)/2, 
        0
    );
}

void SurfaceWindow::wheelEvent(QWheelEvent *event) {
    zoom *= (event->angleDelta().y() > 0) ? 1.1 : 0.9; // Увеличиваем или уменьшаем масштаб
    updateProjection(); // Обновляем проекцию
    update(); // Перерисовываем окно
}

void SurfaceWindow::mousePressEvent(QMouseEvent *event) {
    lastMousePos = event->pos();
}

void SurfaceWindow::mouseMoveEvent(QMouseEvent *event) {
    QPoint delta = event->pos() - lastMousePos;
    if (event->buttons() & Qt::LeftButton) {
        rotationY += delta.x() * 0.5;
        rotationX += delta.y() * 0.5;
        update();
    }
    lastMousePos = event->pos();
}

void SurfaceWindow::drawCoordinateSystem(QPainter &painter) {
    QPen axisPen(Qt::black, 2, Qt::SolidLine);
    painter.setPen(axisPen);

    // Базовые точки осей, соответствующие границам графика
    //QVector3D origin(0, 0, 0);  // Начало координат
    QVector3D x_end(b, 0, 0);   // Конец оси X
    QVector3D y_end(0, d, 0);   // Конец оси Y
    QVector3D z_end(0, 0, m_maxZ);   // Конец оси Z
    QVector3D x_start(a, 0, 0);   // Начало оси X
    QVector3D y_start(0, c, 0);   // Начало оси Y
    QVector3D z_start(0, 0, m_minZ);   // Начало оси Z


    // Проекция базовых точек
    //QVector3D proj_origin = project(origin);
    QVector3D proj_x = project(x_end);
    QVector3D proj_y = project(y_end);
    QVector3D proj_z = project(z_end);
    
    QVector3D proj_x_start = project(x_start);
    QVector3D proj_y_start = project(y_start);
    QVector3D proj_z_start = project(z_start);

    // Рисуем оси
    drawAxis(painter, proj_x_start, proj_x, "X");
    drawAxis(painter, proj_y_start, proj_y, "Y");
    drawAxis(painter, proj_z_start, proj_z, "Z");

    drawCube(painter);
}

void SurfaceWindow::drawAxis(QPainter &painter, const QVector3D &start, const QVector3D &end, const QString &label) {
    // Линия оси
    painter.drawLine(QPointF(start.x(), start.y()), QPointF(end.x(), end.y()));

    // Стрелка
    QPointF arrowP1 = QPointF(end.x(), end.y());
    QVector2D dir(end.x() - start.x(), end.y() - start.y());
    dir = dir.normalized() * 10;
    QPointF arrowP2 = arrowP1 + QPointF(-dir.x() - dir.y(), -dir.y() + dir.x());
    QPointF arrowP3 = arrowP1 + QPointF(-dir.x() + dir.y(), -dir.y() - dir.x());

    QPolygonF arrowHead;
    arrowHead << arrowP1 << arrowP2 << arrowP3;
    painter.drawPolygon(arrowHead);

    // Метка оси
    painter.drawText(QPointF(end.x() + 5, end.y() + 5), label);

    // Метки 0 и 1
    QVector3D zero(0, 0, 0); 
    QVector3D proj_zero = project(zero);

    painter.drawText(QPointF(proj_zero.x() - 15, proj_zero.y() + 5), "0");
    
    QVector3D one_x(1, 0, 0); 
    QVector3D one_y(0, 1, 0);   
    QVector3D one_z(0, 0, 1);   

    QVector3D proj_one_x = project(one_x);
    QVector3D proj_one_y= project(one_y);
    QVector3D proj_one_z = project(one_z);

    if (label == "X") {
        painter.drawText(QPointF(proj_one_x.x() - 15, proj_one_x.y() + 5), QString::number(1));
    } else if (label == "Y") {
        painter.drawText(QPointF(proj_one_y.x() - 15, proj_one_y.y() + 5), QString::number(1));
    } else if (label == "Z") {
        painter.drawText(QPointF(proj_one_z.x() - 15, proj_one_z.y() + 5), QString::number(1));
    }
}

// Функция для рисования куба
void SurfaceWindow::drawCube(QPainter &painter) {
    QVector<QVector3D> cubePoints = {
        QVector3D(a, c, m_minZ), QVector3D(b, c, m_minZ),
        QVector3D(b, d, m_minZ), QVector3D(a, d, m_minZ),
        QVector3D(a, c, m_maxZ), QVector3D(b, c, m_maxZ),
        QVector3D(b, d, m_maxZ), QVector3D(a, d, m_maxZ)
    };

    QVector<QPair<int, int>> cubeEdges = {
        {0,1}, {1,2}, {2,3}, {3,0}, // нижняя грань
        {4,5}, {5,6}, {6,7}, {7,4}, // верхняя грань
        {0,4}, {1,5}, {2,6}, {3,7}  // вертикальные ребра
    };

    QPen cubePen(Qt::gray, 1, Qt::DashLine);
    painter.setPen(cubePen);

    for (auto &edge : cubeEdges) {
        QVector3D p1 = project(cubePoints[edge.first]);
        QVector3D p2 = project(cubePoints[edge.second]);
        painter.drawLine(QPointF(p1.x(), p1.y()),
                         QPointF(p2.x(), p2.y()));
    }
}

void SurfaceWindow::closeEvent(QCloseEvent *event) 
{   
    if (threads_working[0] == true) {
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        event->ignore();  
    } else {
        pthread_mutex_lock(&p_mutex);
        *threads_quiting = true;
        pthread_cond_broadcast(&p_cond);
        pthread_mutex_unlock(&p_mutex);
        event->accept(); 
    }
}
SurfaceWindow::~SurfaceWindow()
{
    for (int i = 0; i < p; i++)
    {
        pthread_join(aA[i].tid, nullptr);
    }
    delete[] aA;
    delete threads_working;
    delete threads_quiting;
}

void SurfaceWindow::clearApproximationData()
{
    approxData.clear();
}

void SurfaceWindow::ApproximationData::clear() 
{   
    if (A) {delete [] A; A = nullptr;}
    if (I) {delete [] I; I = nullptr;}
    if (B) {delete [] B; B = nullptr;}
    if (x) {delete [] x; x = nullptr;}
    if (r) {delete [] r; r = nullptr;}
    if (u) {delete [] u; u = nullptr;}
    if (v) {delete [] v; v = nullptr;}
    calc_status = UNDEF;
}

void SurfaceWindow::ApproximationData::allocate(int nx, int ny, int /*p*/) 
{
    
    N = (nx + 1)*(ny + 1);
    len_msr =  N + 1 + get_len_msr (nx, ny);

    calc_status = UNDEF;
    A = new double[len_msr];
    I = new int[len_msr];
    B = new double[N];
    x = new double[N];
    r = new double[N];
    u = new double[N];
    v = new double[N];
}

void SurfaceWindow::ApproximateFunction()
{
    if(approxData.calc_status == CALCULATED) {return;}
    if(*threads_working == true){return;}
    if(!approxData.A) approxData.allocate(nx, ny, p);

    int k;

    for (k = 0; k < p; k++)
    {
        aA[k].A = approxData.A;
        aA[k].I = approxData.I;
        aA[k].B = approxData.B;
        aA[k].r = approxData.r;
        aA[k].u = approxData.u;
        aA[k].v = approxData.v;
        aA[k].x = approxData.x;
        aA[k].a = a;
        aA[k].b = b;
        aA[k].c = c;
        aA[k].d = d;
        aA[k].func_id = func_id;
        aA[k].maxit = max_it;
        aA[k].p = p;
        aA[k].k = k;
        aA[k].eps = eps;
        aA[k].nx = nx;
        aA[k].ny = ny;
        aA[k].N = (nx + 1)*(ny + 1);
        aA[k].working = threads_working;
        aA[k].quiting_app = threads_quiting;
        aA[k].len_msr = approxData.len_msr;
        aA[k].p_mutex = &p_mutex;
        aA[k].p_cond = &p_cond;
        aA[k].norm = norm;
        aA[k].point = point;
    }

    pthread_mutex_lock(&p_mutex);
    *threads_working = true;
    pthread_cond_broadcast(&p_cond);
    pthread_mutex_unlock(&p_mutex);
}

SurfaceWindow::SurfaceWindow(QWidget *parent) : QWidget(parent) {
    setFocusPolicy(Qt::StrongFocus);
    updateProjection();
    a = -1, b = 1, c = -1, d = 1, eps = 1e-10;
    nx = 5, ny = 5, func_id = 7, max_it = 100, p=1;
    first = true;
    point = 0;
    norm = 0;
    threads_working = new bool(false);
    threads_quiting = new bool(false);
    currentFunc = FUNC;
    max = 1;
    strochka_printed = false;


    change_func();
}

int SurfaceWindow::parse_command_line(int argc, char* argv[])
{

    if (argc != 13 || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1 || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1|| sscanf(argv[7], "%d", &mx) != 1 || sscanf(argv[8], "%d", &my) != 1 || sscanf(argv[9], "%d", &func_id) != 1 || sscanf(argv[10], "%lf", &eps) != 1 || sscanf(argv[11], "%d", &max_it) != 1 || sscanf(argv[12], "%d", &p) != 1) 
    {
        printf("Usage ./a.out a b c d nx ny mx my func_id epsilon max_iterations p\n");
        return 1;   
    }
    if(func_id < 0 || nx < 0|| ny < 0 || max_it < 0 || p < 0 || func_id > 7)
    {
        printf("Usage ./a.out a b c d nx ny mx my func_id epsilon max_iterations p\n");
        return 1;   
    }
    
    aA = new Args[p];
    
    approxData.allocate(nx, ny, p);

    int k;
    *threads_working = false;
    for (k = 0; k < p; k++)
    {
        aA[k].A = approxData.A;
        aA[k].I = approxData.I;
        aA[k].B = approxData.B;
        aA[k].r = approxData.r;
        aA[k].u = approxData.u;
        aA[k].v = approxData.v;
        aA[k].x = approxData.x;
        aA[k].a = a;
        aA[k].b = b;
        aA[k].c = c;
        aA[k].d = d;
        aA[k].func_id = func_id;
        aA[k].maxit = max_it;
        aA[k].p = p;
        aA[k].k = k;
        aA[k].eps = eps;
        aA[k].nx = nx;
        aA[k].ny = ny;
        aA[k].N = (nx + 1)*(ny + 1);
        aA[k].working = threads_working;
        aA[k].quiting_app = threads_quiting;
        aA[k].len_msr = approxData.len_msr;
        aA[k].p_mutex = &p_mutex;
        aA[k].p_cond = &p_cond;
        aA[k].norm = norm;
        aA[k].point = point;
    }

    for (k = 0; k < p; k++)
    {
        if (pthread_create(&aA[k].tid, nullptr, thread_func, aA + k))
        {
            std::cerr << "Error creating thread " << k << std::endl;
            clearApproximationData();
            return 1;
        }
    }
    return 0;
}