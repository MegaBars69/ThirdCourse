#include "surfacewindow.h"
#include <cmath>
#include <algorithm>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QPainterPath>
#include <QMessageBox>
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


int SurfaceWindow::parse_command_line(int argc, char* argv[])
{

    if (argc != 13 || sscanf(argv[1], "%lf", &a) != 1 || sscanf(argv[2], "%lf", &b) != 1 || sscanf(argv[3], "%lf", &c) != 1 || sscanf(argv[4], "%lf", &d) != 1 || sscanf(argv[5], "%d", &nx) != 1 || sscanf(argv[6], "%d", &ny) != 1|| sscanf(argv[7], "%d", &mx) != 1 || sscanf(argv[8], "%d", &my) != 1 || sscanf(argv[9], "%d", &func_id) != 1 || sscanf(argv[10], "%lf", &eps) != 1 || sscanf(argv[11], "%d", &max_it) != 1 || sscanf(argv[12], "%d", &p) != 1) 
    {
        printf("Usage1 ./a.out a b c d nx ny mx my func_id epsilon max_iterations p\n");
        return 1;   
    }
    if(func_id < 0 || nx < 0|| ny < 0 || max_it < 0 || p < 0 || func_id > 7)
    {
        printf("Usage1 ./a.out a b c d nx ny mx my func_id epsilon max_iterations p\n");
        return 1;   
    }
    return 0;
}

void SurfaceWindow::change_func()
{
    if (!first)
    {
        func_id = (func_id + 1) % 8;
    }
    else
    {
        first = false;
    }
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
    if (currentFunc == FUNC)
        currentFunc = APPROX;
    else if (currentFunc == APPROX)
        currentFunc = ERRORS;
    else
        currentFunc = FUNC;

    update();
}


void SurfaceWindow::zoom_in()
{
    double center = (a + b) / 2.0;
    double half_length = (b - a) / 4.0;
    a = center - half_length;
    b = center + half_length;
    
    center = (c + d) / 2.0;
    half_length = (d - c) / 4.0;
    c = center - half_length;
    d = center + half_length;
    clearApproximationData();
    update();
}

void SurfaceWindow::zoom_out()
{
    double center = (a + b) / 2.0;
    double half_length = (b - a);
    a = center - half_length;
    b = center + half_length;

    center = (c + d) / 2.0;
    half_length = (d - c);
    c = center - half_length;
    d = center + half_length;
    clearApproximationData();
    update();
}

void SurfaceWindow::increase_points()
{
    nx *= 2;
    ny *= 2;
    clearApproximationData();
    update();
}

void SurfaceWindow::decrease_points()
{
    if (nx > 3 && ny > 3)
    {
        nx /= 2;
        ny /= 2;
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
    point--;
    clearApproximationData();
    update();
}

void SurfaceWindow::point_up()
{
    point++;
    clearApproximationData();
    update();
}

void SurfaceWindow::keyPressEvent(QKeyEvent *event) {
    // Проверяем статус вычислений
    //pthread_mutex_lock(&p_mutex);    
    if (approxData.calc_status == CALCULATING) 
    {
        printf("В данный момент производятся вычисления. Пожалуйста, подождите.\n");
        QMessageBox::information(this, 
                                "Вычисление", 
                                "В данный момент производятся вычисления. Пожалуйста, подождите.");
        return;  // Прерываем обработку нажатия
    }
    else if (approxData.calc_status == CALCULATED) {
        std::cout<<"Calculated'\n";
    }
    else{
        std::cout<<"Undef'\n";
    }
    //pthread_mutex_unlock(&p_mutex);


    // Если вычислений нет - обрабатываем клавиши как обычно
    //pthread_mutex_lock(&p_mutex);
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
        default: QWidget::keyPressEvent(event);
    }
    /*pthread_mutex_unlock(&p_mutex);
    pthread_cond_broadcast(&p_cond);*/
}

void SurfaceWindow::reset_view() {
    zoom = 1.33;
    rotationX = -68.5;
    rotationY = -0.5;
    update();
}

SurfaceWindow::SurfaceWindow(QWidget *parent) : QWidget(parent) {
    setFocusPolicy(Qt::StrongFocus);
    updateProjection();
    a = -1, b = 1, c = -1, d = 1, eps = 1e-10;
    nx = 5, ny = 5, func_id = 7, max_it = 100, p=1;
    first = true;
    point = 0;
    currentFunc = FUNC;
    change_func();
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

void SurfaceWindow::calculateSurface() 
{
    vertices.clear();
    triangles.clear();

    double xi, yj, z = 0;

    double hx_loc = (b - a) / mx;
    double hy_loc = (d - c) / my;
    double hx = (b - a) / nx;
    double hy = (d - c) / ny;
    // Инициализация min/max
    m_minZ = std::numeric_limits<double>::max();
    m_maxZ = std::numeric_limits<double>::lowest();
    if(currentFunc != FUNC && approxData.calc_status != CALCULATED)
    {
        approxData.calc_status = CALCULATING;
        ApproximateFunction();

        approxData.calc_status = CALCULATED;
    }

    // Generate vertices
    for (int i = 0; i <= mx; ++i) {
        for (int j = 0; j <= my; ++j) {
            xi = a + i * hx_loc;
            yj = c + j * hy_loc;
            if(currentFunc == FUNC)
            {
                z = f(xi,yj);
            }
            else if(currentFunc == APPROX && approxData.calc_status == CALCULATED)
            {
                z = Pf(approxData.x, xi, yj, a, c, hx, hy, nx, ny);
            }
            else if(currentFunc == ERRORS && approxData.calc_status == CALCULATED)
            {
                z = fabs(f(xi,yj) - Pf(approxData.x, xi, yj, a, c, hx, hy, nx, ny));
            }
                
            // Обновляем min/max
            if (z < m_minZ) m_minZ = z;
            if (z > m_maxZ) m_maxZ = z;
            norm = (fabs(m_maxZ) > fabs(m_minZ) ? fabs(m_maxZ) : fabs(m_minZ));

            vertices.append(QVector3D(xi, yj, z));
        }
    }

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
    //update();
}
void SurfaceWindow::paintEvent(QPaintEvent *) {
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
    // Sort triangles by depth
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
        
        // Нормализованная высота для цвета
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
    sprintf(strochka, "%s%e", prefix, norm);
    painter.drawText(0, 130, strochka);
    
    prefix = "(nx, ny) = ";
    sprintf(strochka, "%s%d", prefix, nx);
    sprintf(strochka, "%s%d", prefix, ny);
    painter.drawText(0, 95, strochka);
    
    prefix = "(mx, my) = ";
    sprintf(strochka, "%s%d", prefix, mx);
    sprintf(strochka, "%s%d", prefix, my);
    painter.drawText(0, 115, strochka);

    prefix = "p = ";
    sprintf(strochka, "%s%d", prefix, p);
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
    delete[] strochka;
}

QVector3D SurfaceWindow::project(const QVector3D &point) const {
    QMatrix4x4 view;
    view.translate(0, 0, -5 * zoom);
    view.rotate(rotationX, 1, 0, 0);
    view.rotate(rotationY, 0, 1, 0);
    
    QVector3D proj = projection * view * point;
    proj.setX((proj.x() + 1) * width() / 2);
    proj.setY((1 - proj.y()) * height() / 2);
    return proj;
}

void SurfaceWindow::updateProjection()
{
    projection.setToIdentity();
    
    // Рассчитываем границы с учетом Z
    QSize viewport = size();
    double aspect = static_cast<double>(viewport.width()) / viewport.height();
    
    // Ортографическая проекция с автоматическим масштабированием
    QVector3D sceneSize(b - a, d - c, m_maxZ - m_minZ);
    double maxSceneSize = qMax(sceneSize.x(), qMax(sceneSize.y(), sceneSize.z()));
    
    // Применяем масштабирование (zoom)
    maxSceneSize *= zoom;
    
    projection.ortho(
        -maxSceneSize * aspect/2, maxSceneSize * aspect/2, // left/right
        -maxSceneSize/2,         maxSceneSize/2,          // bottom/top
        -maxSceneSize*2,          maxSceneSize*2           // near/far
    );
    
    // Центрируем сцену
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

    // Рисуем куб (опционально)
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

void SurfaceWindow::clearApproximationData()
{
    approxData.clear();
}

void SurfaceWindow::ApproximationData::clear() 
{   
    /*if(calc_status != CALCULATING)
    {*/
        if (A) {delete [] A; A = nullptr;}
        if (I) {delete [] I; I = nullptr;}
        if (B) {delete [] B; B = nullptr;}
        if (x) {delete [] x; x = nullptr;}
        if (r) {delete [] r; r = nullptr;}
        if (u) {delete [] u; u = nullptr;}
        if (v) {delete [] v; v = nullptr;}
        if (aA) {delete[] aA; aA = nullptr;}
        free_reduce_sum();
        calc_status = UNDEF;
    //}
}

void SurfaceWindow::ApproximationData::allocate(int nx, int ny, int p) 
{
    //clear();
   /*if(calc_status == UNDEF)
    {*/
        N = (nx + 1)*(ny + 1);
        len_msr =  N + 1 + get_len_msr (nx, ny);

        calc_status = UNDEF;
        aA = new Args[p];
        A = new double[len_msr];
        I = new int[len_msr];
        B = new double[N];
        x = new double[N];
        r = new double[N];
        u = new double[N];
        v = new double[N];
        init_reduce_sum (p);
    //}
}

void SurfaceWindow::ApproximateFunction()
{
    if(approxData.calc_status == CALCULATED) {return;}
    if(!approxData.A) approxData.allocate(nx, ny, p);

    Args* aA = approxData.aA;
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
        aA[k].len_msr = approxData.len_msr;
    }

    //Запуск потоков.
    for (k = 0; k < p; k++)
    {
        if (pthread_create(&aA[k].tid, nullptr, thread_func, aA+k))
        {
            std::cerr << "Error creating thread " << k <<std::endl;
            clearApproximationData();
            return;
        }
    }
    
    //Прибиваем потоки молотком.
    for (k = 0; k < p; k++)
    {
        pthread_join(aA[k].tid, nullptr);
    }

    its = aA->its;
    r1 = aA->r1;
    r2 = aA->r2;
    r3 = aA->r3;
    r4 = aA->r4;
    t1 = aA->t1;
    t2 = aA->t2;
    printf (
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
        It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        "./a.out", 5, r1, r2, r3, r4, t1, t2, its, eps, func_id, nx, ny, p);
    
}