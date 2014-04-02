QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QMAKE_CXXFLAGS += -std=c++11 -pedantic

TEMPLATE = lib
CONFIG += staticlib create_prl c++11 link_prl
DEPENDPATH += .
INCLUDEPATH += ..

# Input
HEADERS += calculations.h optimization_task.h processing.h

SOURCES += calculations.cpp optimization_task.cpp processing.cpp

LIBS += -L/usr/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui \
	-L../cpp_utils -lcpp_utils \
	-L../qt_utils -lqt_utils \

