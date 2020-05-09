TEMPLATE = app

CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS_DEBUG += -O3 -fopenmp -larmadillo

SOURCES += \
    ../Code_Files/Pro5_Functions.cpp \
    ../Code_Files/Pro5_Tests.cpp \
    ../Code_Files/Pro5_main.cpp


INCLUDEPATH += C:\Qt\armadillo-9.700.2\include
DEPENDPATH += C:\Qt\armadillo-9.700.2\include


LIBS += \
    -LC:\Qt\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT \
    -fopenmp

HEADERS += \
    ../Code_Files/Pro5_Functions.h
