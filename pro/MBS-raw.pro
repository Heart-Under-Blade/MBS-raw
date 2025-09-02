TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle

DESTDIR = ../bin

VERSION = 3.0.0

QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS += -std=gnu++11
QMAKE_CXXFLAGS += -march=corei7 -msse4.2

CONFIG(release, debug|release): {
    TARGET = mbs
}
CONFIG(debug, debug|release): {
    DEFINES += _DEBUG
    TARGET = mbs_d
}

SRC = ../src/main

INCLUDEPATH += \
    $$SRC \
    $$SRC/math \
    ../src/bigint \
    $$SRC/handler \
    $$SRC/service \
    $$SRC/geometry \
    $$SRC/geometry/sse \
    $$SRC/splitting

SOURCES += \
    $$files($$SRC/*.cpp, true) \
    $$files(../src/bigint/*.cc, true)

HEADERS += \
    $$files($$SRC/*.h, true) \
    $$files(../src/bigint/*.hh, true) \

#message($$SOURCES)
