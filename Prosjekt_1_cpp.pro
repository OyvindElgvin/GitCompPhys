TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


INCLUDEPATH += /usr/local/Cellar/armadillo/9.600.6/include/
LIBS += -L/usr/local/Cellar/armadillo/9.600.6/lib/ -larmadillo


SOURCES += \
        main.cpp
