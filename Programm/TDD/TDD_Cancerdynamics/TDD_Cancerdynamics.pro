#-------------------------------------------------
#
# Project created by QtCreator 2014-09-17T10:54:09
#
#-------------------------------------------------

QT       += testlib

QMAKE_CXXFLAGS += -std=c++0x
CONFIG += c++0x
QMAKE_MAC_SDK = macosx10.9

QT       -= gui

TARGET = tst_tdd_cancerdynamicstest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += tst_tdd_cancerdynamicstest.cpp \
    ../../Simulator/Cancerdynamics/eventclass.cpp \
    ../../Simulator/Cancerdynamics/graphclass.cpp \
    ../../Simulator/Cancerdynamics/randomdice.cpp \
    ../../Simulator/Cancerdynamics/traitclass.cpp \
    ../../Simulator/Cancerdynamics/traiteventmanager.cpp \
    Testing_Utilities/Utilities.cpp \
    ../../Simulator/Cancerdynamics/FileStreamer/cfilestreamer.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../../Simulator/Cancerdynamics/eventclass.h \
    ../../Simulator/Cancerdynamics/graphclass.h \
    ../../Simulator/Cancerdynamics/randomdice.h \
    ../../Simulator/Cancerdynamics/traitclass.h \
    ../../Simulator/Cancerdynamics/traiteventmanager.h \
    ../../Simulator/Cancerdynamics/FileStreamer/cfilestreamer.h \
    ../../Simulator/Cancerdynamics/FileStreamer/ifilestreamer.h \
    ../../Simulator/Cancerdynamics/FileStreamer/operators.h \
    Testing_Utilities/Utilities.h
