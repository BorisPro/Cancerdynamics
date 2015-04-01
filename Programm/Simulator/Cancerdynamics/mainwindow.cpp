#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "filestreaming.h"
#include <QFileDialog>
#include <thread>
#include <cmath>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    Plot(new PlotWindow)
{
    ui->setupUi(this);
//    ui->pushButton_plot->setEnabled(true);
    /// connect the the signal clicked() from Button"plot" to the SLOT
    /// show() from the Plot Window
    connect(ui->pushButton_plot, SIGNAL(clicked()),Plot,SLOT(show()));
    connect(Plot,SIGNAL(graphDrawn()),this,SLOT(enablePlotButton()));
}

MainWindow::~MainWindow()
{
    delete Plot;
    delete ui;
}

void MainWindow::createBlankInstance()
{

}

void MainWindow::addTreePopulationProperties()
{
//    addRootItem("Population Properties");
//    addChildItem(ui->treeWidget_parameters->itemAt(0,0) ,"Size: " + QString::number(TraitClass::Size));
//    addChildItem(ui->treeWidget_parameters->itemAt(0,0) ,"K:    " + QString::number(TraitClass::K));
//    addChildItem(ui->treeWidget_parameters->itemAt(0,0) ,"Mutation (K ingored): " + QString::number(TraitClass::Mutation * (TraitClass::K * sqrt(TraitClass::K))));
//    addChildItem(ui->treeWidget_parameters->itemAt(0,0) ,"Competition (K ingored)");
//    for(int i = 0; i < TraitClass::Size; i++){
//        QString row("");
//        for(int j = 0; j < TraitClass::Size; j++)
//            row +=QString::number(TraitClass::CompDeathRate[i][j]*TraitClass::K) + " ";
//        addChildItem(ui->treeWidget_parameters->itemAt(0,0)->child(3), row);
//    }
}

void MainWindow::addTreeTraitProperties(PopulationManager Parameters)
{
//    addRootItem("Members  (K ingored)");
//    addRootItem("Birth Rates");
//    addRootItem("Death Rates");
//    for(int i = 0; i < TraitClass::Size; i++){
//        addChildItem(ui->treeWidget_parameters->topLevelItem(1)
//                     ,QString::number(i+1) + ". | " + QString::number(Parameters.getKMembers(i)));
//        addChildItem(ui->treeWidget_parameters->topLevelItem(2)
//                     ,QString::number(i+1) + ". | " + QString::number(Parameters.Trait[i].BirthRate));
//        addChildItem(ui->treeWidget_parameters->topLevelItem(3)
//                     ,QString::number(i+1) + ". | " + QString::number(Parameters.Trait[i].DeathRate));
//    }
}

void MainWindow::displayTraitData(PopulationManager Parameters)
{
//    ui->treeWidget_parameters->clear();
//    ui->treeWidget_parameters->setEnabled(isFound);
//    ui->groupBox_Plot->setEnabled(isFound);
//    if(isFound){
//        addTreePopulationProperties();
//        addTreeTraitProperties(Parameters);
//    }
//    else
//        addRootItem("file not found!");

//    ui->tableWidget_Fitness->clear();
//    ui->tableWidget_Fitness->setColumnCount(TraitClass::Size);
//    ui->tableWidget_Fitness->setRowCount(TraitClass::Size);
//    // set the default size, here i've set it to 20px by 20x
//    ui->tableWidget_Fitness->horizontalHeader()->setDefaultSectionSize(40);
//    ui->tableWidget_Fitness->verticalHeader()->setDefaultSectionSize(40);
//    for(int i = 0; i < TraitClass::Size; ++i){
//        for(int j = 0; j < TraitClass::Size; ++j){
//            if(std::abs(i-j) <= 1)
//                ui->tableWidget_Fitness->setItem(i,j, new QTableWidgetItem(QString::number(TraitClass::Fitness[i][j])));
//            else
//                ui->tableWidget_Fitness->setItem(i,j, new QTableWidgetItem("-"));
//        }
//    }
//    for(int i = 0; i < TraitClass::Size; ++i){
//        for(int j = 0; j < TraitClass::Size; ++j){
//            if(std::abs(i-j) > 1){
//                ui->tableWidget_Fitness->item(j,i)->setBackgroundColor(QColor(Qt::black));
//                continue;
//            }
//            else if(TraitClass::Fitness[i][j]*TraitClass::Fitness[j][i] >= 0 && i != j)
//                ui->tableWidget_Fitness->item(i,j)->setBackgroundColor(QColor(Qt::red));
//            else if(std::max(TraitClass::Fitness[j][i]/Parameters.Trait[j].BirthRate,0.) > 0.5)
//                ui->tableWidget_Fitness->item(j,i)->setBackgroundColor(QColor(Qt::green));
//        }
//    }
}

void MainWindow::iterateMembers(QString StepName, int size)
{
    for(int i = 1; i < size; i++)
        if(StepName == "enter " + QString::number(i) + ". Members:")
            ui->label_InputName->setText("enter " + QString::number(i+1) + ". Members:");
    if(StepName == "enter " + QString::number(size) + ". Members:")
        ui->label_InputName->setText("enter 1. Births:");
}

void MainWindow::iterateBirths(QString StepName, int size)
{
    for(int i = 1; i < size; i++)
        if(StepName == "enter " + QString::number(i) + ". Births:")
            ui->label_InputName->setText("enter " + QString::number(i+1) + ". Births:");
    if(StepName == "enter " + QString::number(size) + ". Births:")
        ui->label_InputName->setText("enter 1. Deaths:");
}

void MainWindow::iterateDeaths(QString StepName, int size)
{
    for(int i = 1; i < size; i++)
        if(StepName == "enter " + QString::number(i) + ". Deaths:")
            ui->label_InputName->setText("enter " + QString::number(i+1) + ". Deaths:");
    if(StepName == "enter " + QString::number(size) + ". Deaths:"){
        ui->label_InputName->setText("file created!");
        ui->groupBox_Input->setEnabled(false);
        ui->groupBox_FileStream->setEnabled(true);
    }
}

void MainWindow::iterateInputParameter()
{
//    QString StepName = ui->label_InputName->text();
//    QString value = ui->lineEdit_Input->text();
//    addRootItem(value);
//    FileStreaming::appendToFileEndl(ui->lineEdit_FileName->text(), value);
//    int size = (ui->treeWidget_parameters->topLevelItem(0)->text(0)).toInt();

//    iterateTraitProperties(StepName, size);
//    iterateMembers(StepName, size);
//    iterateBirths(StepName, size);
//    iterateDeaths(StepName, size);
}

bool MainWindow::iterateTraitProperties(QString StepName, int size)
{
    if(StepName == "enter Size:")
        ui->label_InputName->setText("enter K:");
    else if(StepName == "enter K:")
        ui->label_InputName->setText("enter Mutationprobability:");
    else if(StepName == "enter Mutationprobability:")
        ui->label_InputName->setText("enter 1. line of competition:");
    for(int i = 1; i < size; i++){
        if(StepName == "enter " + QString::number(i) + ". line of competition:")
            ui->label_InputName->setText("enter " + QString::number(i+1) + ". line of competition:");
    }
    if(StepName == "enter " + QString::number(size) + ". line of competition:")
        ui->label_InputName->setText("enter 1. Members:");
    return true;
}

void MainWindow::readFileContent(QString Filename){
    try{
        PopulationManager Parameters;
        Parameters.initWithFile(Filename.toStdString());
        displayTraitData(Parameters);
    }
    catch(std::string exception){
        QMessageBox msgBox;
        msgBox.setText(QString::fromStdString(exception));
        msgBox.exec();
    }
}

void MainWindow::addRootItem(QString value)
{
    QTreeWidgetItem* root = new QTreeWidgetItem(ui->treeWidget_parameters);
    root->setText(0, value);
    ui->treeWidget_parameters->addTopLevelItem(root);
}

void MainWindow::addChildItem(QTreeWidgetItem *parent, QString value)
{
    QTreeWidgetItem* root = new QTreeWidgetItem();
    root->setText(0,value);
    parent->addChild(root);
}

// Slots ----------------------------------------------

void MainWindow::on_pushButton_plot_clicked()
{
    int iterations = ui->lineEdit_maxIterations->text().toInt();
    QString FName = ui->lineEdit_FileName->text();
//    if(!FName.contains("/Parameters/"))
//        FName = "../Paramters/"+FName;
    try{
        Plot->rePlot(iterations, FName);
    }
    catch(std::string exception){
        QMessageBox msgBox;
        msgBox.setText(QString::fromStdString(exception));
        msgBox.exec();
    }
    ui->pushButton_plot->setEnabled(false);
}

void MainWindow::on_pushButton_load_File_clicked()
{
//    QString Filename = ui->lineEdit_FileName->text();
//    ui->treeWidget_parameters->setHeaderLabel(Filename);
//    readFileContent(Filename);
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                     "",
                                                     tr("Files (*)"));
    ui->lineEdit_FileName->setText(fileName);
//    ui->pushButton_plot->setEnabled(false);
//    if(fileName != ""){
//        ui->pushButton_plot->setEnabled(true);
//    }
}

void MainWindow::on_pushButton_create_File_clicked()
{
//    ui->treeWidget_parameters->clear();
//    ui->treeWidget_parameters->setEnabled(true);
//    ui->groupBox_Input->setEnabled(true);
//    ui->groupBox_FileStream->setEnabled(false);
//    ui->label_InputName->setText("enter Size:");
//    QString Filename = ui->lineEdit_FileName->text();
//    ui->treeWidget_parameters->setHeaderLabel(Filename);
}

void MainWindow::on_lineEdit_Input_returnPressed()
{
    iterateInputParameter();
    ui->lineEdit_Input->clear();
}

void MainWindow::enablePlotButton()
{
    ui->pushButton_plot->setEnabled(true);
}
