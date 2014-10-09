#-------------------------------------------------
#
# Project created by QtCreator 2014-09-16T10:53:57
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Cancerdynamics
TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++0x
CONFIG += c++0x
QMAKE_MAC_SDK = macosx10.9

SOURCES += main.cpp\
        mainwindow.cpp \
    eventclass.cpp \
    graphclass.cpp \
    plotwindow.cpp \
    randomdice.cpp \
    traitclass.cpp \
    traiteventmanager.cpp \
    qcustomplot.cpp \
    FileStreamer/cfilestreamer.cpp

HEADERS  += mainwindow.h \
    eventclass.h \
    graphclass.h \
    plotwindow.h \
    randomdice.h \
    traitclass.h \
    traiteventmanager.h \
    qcustomplot.h \
    FileStreamer/cfilestreamer.h \
    FileStreamer/ifilestreamer.h \
    FileStreamer/operators.h

FORMS    += mainwindow.ui \
    plotwindow.ui
