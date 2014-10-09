#include "plotwindow.h"
#include "ui_plotwindow.h"
#include <cmath>
#include <QThread>

bool isNear(double Value, double Expected, double K){
    double diff = Value - Expected;
    if( diff < 10./K && diff > -10./K)
        return true;
    return false;
}

/** \class PlotRenderer
 * \brief Simple QObject for data rendering
 *
 * This class takes account of the heavy work
 * which is being done by GraphClass::generateEvolution().
 *
 * Since generateEvolution() can take a huge amount of time,
 * one can simply move this class to another QThread
 *
 * \sa PlotWindow
*/

// --------------------------------------

PlotWindow::PlotWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PlotWindow)
{
    ui->setupUi(this);
    ui->customPlot->xAxis->setLabel("Time");
    ui->customPlot->yAxis->setLabel("Cells");

    PlotRenderer * plotter = new PlotRenderer;

    QThread * thread = new QThread(this);
    connect(this,SIGNAL(destroyed()), thread, SLOT(quit()));

    plotter->moveToThread(thread);

    connect(this,   SIGNAL(request_plot(int,QString)),
            plotter,SLOT(renderData(int,QString)),
            Qt::QueuedConnection);
    connect(plotter,SIGNAL(renderedData(const GraphClass*,int)),
            this, SLOT(drawGraph(const GraphClass*,int)),
            Qt::QueuedConnection);
    thread->start();
}

void PlotWindow::addGraphWithName(QString GName, QPen graphPen)
{
    ui->customPlot->addGraph();
    ui->customPlot->graph()->setName(GName);
    ui->customPlot->graph()->setLineStyle(QCPGraph::lsStepLeft);
    ui->customPlot->graph()->setPen(graphPen);
}

void PlotWindow::createGraphs()
{
    QVector<QPen> graphPen(numberOfGraphs);
    for(size_t currentTrait = 0; currentTrait < numberOfGraphs; currentTrait++){
        graphPen[currentTrait].setColor(QColor(rand()%245+10, rand()%245+10, rand()%245+10));
        addGraphWithName(QString::number(currentTrait+1) + ". Trait", graphPen[currentTrait]);
    }
    ui->customPlot->legend->setVisible(true);
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

void PlotWindow::fillCreatedGraphs(const GraphClass * Graph)
{
    ui->customPlot->yAxis->setRange(0,Graph->getMaxMembers());
    ui->customPlot->xAxis->setRange(0,Graph->getMaxTime());
    for(size_t i = 0; i < numberOfGraphs; ++i){
        ui->customPlot->graph(i)->addData(Graph->getTimesOf(i),Graph->getTraitHistOf(i));
        ui->customPlot->graph(i)->setData(Graph->getTimesOf(i),Graph->getTraitHistOf(i));
    }
    ui->customPlot->replot();
    delete Graph;
}

void PlotWindow::rePlot(int iterations, QString FName)
{
    ui->lineEdit_ImageName->setText(FName);
    emit request_plot(iterations, FName);   // wird empfangen von "renderData(int,QString)"
}

/** \brief Renders the data given by \a FName
 * \param iterations max number of iterations
 * \param FName file name of the data
 * \param RangeCheck whether to check the range of the data
 */
void PlotRenderer::renderData(int iterations, QString FName){
    //! We use GraphClass to create the data
    try{
        GraphClass * graph = new GraphClass(FName);
        int it_made = graph->generateEvolution(iterations);
        emit renderedData(graph, it_made); // wird empfangen von "drawGraph(const GraphClass*,int)"
    }
    catch(std::string exception){
        QMessageBox msgBox;
        msgBox.setText(QString::fromStdString(exception));
        msgBox.exec();
    }
}

void PlotWindow::drawGraph(const GraphClass * Graph, int iterations){
    try{
        numberOfGraphs = Graph->getNumberOfGraphs();
        ui->customPlot->clearGraphs();
        ui->label_iterations->setText("Iterations: " + QString::number(iterations));
        createGraphs();
        fillCreatedGraphs(Graph);
        emit graphDrawn(); // wird vom MainWindow empfangen und gibt den Plot Button frei.
    }
    catch(std::string exception){
        QMessageBox msgBox;
        msgBox.setText(QString::fromStdString(exception));
        msgBox.exec();
    }
}

void PlotWindow::on_pushButton_saveImage_clicked()
{
    ui->customPlot->savePdf(/*"../Parameters/" + */ui->lineEdit_ImageName->text()+ ".pdf");
    ui->customPlot->savePng(/*"../Parameters/" + */ui->lineEdit_ImageName->text()+ ".png");
}

PlotWindow::~PlotWindow()
{
    delete ui;
}

