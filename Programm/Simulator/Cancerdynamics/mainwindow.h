#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCore>
#include <QtGui>
#include <QMessageBox>
#include "plotwindow.h"
#include "ui_mainwindow.h"
#include <QTextStream>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void readFileContent(QString Filename);
    ~MainWindow();

    // creation of new cancer dynamics instances
    void createBlankInstance();

    // manages the display of data in the tree widget!
    void addTreePopulationProperties();
    void addTreeTraitProperties(PopulationManager Parameters);
    void displayTraitData(/*QString Filename, */PopulationManager Parameters);
    void addRootItem(QString value);
    void addChildItem(QTreeWidgetItem* parent, QString value);

    // manages the creation of new instances.
    void iterateInputParameter();
    bool iterateTraitProperties(QString StepName, int size);
    void iterateMembers(QString StepName, int size);
    void iterateBirths(QString StepName, int size);
    void iterateDeaths(QString StepName, int size);


private slots:
    void on_pushButton_plot_clicked();

    void on_pushButton_load_File_clicked();

    void on_pushButton_create_File_clicked();

    void on_lineEdit_Input_returnPressed();

    void enablePlotButton();    // recieves signal "graphDrawn"
private:
    Ui::MainWindow *ui;
    PlotWindow * Plot;
};

#endif // MAINWINDOW_H
