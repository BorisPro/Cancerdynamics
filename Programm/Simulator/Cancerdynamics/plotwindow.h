#ifndef PLOTWINDOW_H
#define PLOTWINDOW_H

#include <QMessageBox>
#include <QWidget>
#include "qcustomplot.h"
#include "graphclass.h"

namespace Ui {
    class PlotWindow;
}

class PlotRenderer : public QObject {
    Q_OBJECT
signals:
    void renderedData(const GraphClass *, int);

public slots:
    void renderData(int iterations, QString FName);
};

class PlotWindow : public QWidget
{
    Q_OBJECT

signals:
    /// Es wird ein Plot angefertigt, der mit "QString" Parametern generiert wird und "int" Punkt enthalten soll.
    void request_plot(int, QString);
    /// Die Berechnungen sind abgeschlossen und der Graph fertig vorbereitet.
    void graphDrawn();

public:
    explicit PlotWindow(QWidget *parent = 0);
    ~PlotWindow();
    void rePlot(int iterations, QString FName);

private:
    void createGraphs();
    void fillCreatedGraphs(const GraphClass *Graph);
    void addGraphWithName(QString GName, QPen graphPen);


public slots:
    void drawGraph(const GraphClass *g, int iterations);


private slots:
    void on_pushButton_saveImage_clicked();

private:
    size_t numberOfGraphs;
    Ui::PlotWindow *ui;
};

#endif // PLOTWINDOW_H
