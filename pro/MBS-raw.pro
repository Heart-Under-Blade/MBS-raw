TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle

DESTDIR = ../bin

VERSION = 3.0.1

QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2
QMAKE_CXXFLAGS += -Wall -Wextra

LIBS += -lm

CONFIG(release, debug|release): {
    QMAKE_CXXFLAGS += -O2 -flto
    QMAKE_LFLAGS += -flto
    TARGET = mbs
}

CONFIG(debug, debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_d
}

INCLUDEPATH += \
    ../src \
    ../src/math \
    ../src/common \
    ../src/particle \
    ../src/geometry \
    ../src/scattering \
    ../src/handler \
    ../src/tracer \
    ../src/bigint

SOURCES += \
    ../src/*.cpp \
    ../src/math/*.cpp \
    ../src/particle/*.cpp \
    ../src/geometry/*.cpp \
    ../src/handler/*.cpp \
    ../src/tracer/*.cpp \
    ../src/common/*.cpp \
    ../src/scattering/*.cpp \
    ../src/bigint/*.cc

HEADERS += \
    ../src/*.h \
    ../src/math/*.hpp \
    ../src/math/*.h \
    ../src/particle/*.h \
    ../src/geometry/*.h \
    ../src/handler/*.h \
    ../src/tracer/*.h \
    ../src/common/*.h \
    ../src/scattering/*.h \
    ../src/bigint/*.hh
