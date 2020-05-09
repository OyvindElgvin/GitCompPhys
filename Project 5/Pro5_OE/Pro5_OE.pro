TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#CONFIG += -lgomp


INCLUDEPATH += /usr/local/Cellar/armadillo/9.600.6/include/
LIBS += -L/usr/local/Cellar/armadillo/9.600.6/lib/ -larmadillo
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib


#LIBS += -llapack -lblas
#TARGET = OpenMp
#LIBS += -lgomp

QMAKE_CXXFLAGS_DEBUG += -lomp
QMAKE_LFLAGS += -lomp
#QMAKE_CXXFLAGS_DEBUG += -fopenmp
#LIBS += -fopenmp
LIBS += -lomp



SOURCES += \
        ../Code_Files/Pro5_Functions.cpp \
        ../Code_Files/Pro5_Tests.cpp \
        ../Code_Files/Pro5_main.cpp

HEADERS += \
    ../Code_Files/Pro5_Functions.h










