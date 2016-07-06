TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Delmesh.cpp \
    geom_2d.cpp \
    triangle_properties.cpp

HEADERS += \
    Delmesh.h \
    priorityqueue.h \
    mytraits.h \
    geom_2d.h \
    triangle_properties.hpp

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/ -lOpenMeshCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/ -lOpenMeshCore
else:unix: LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/ -lOpenMeshCore

INCLUDEPATH += $$PWD/../../programs/OpenMesh-5.1/build/Build
DEPENDPATH += $$PWD/../../programs/OpenMesh-5.1/build/Build

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/libOpenMeshCore.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/libOpenMeshCore.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/OpenMeshCore.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/OpenMeshCore.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/libOpenMeshCore.a

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/ -lOpenMeshTools
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/ -lOpenMeshTools
else:unix: LIBS += -L$$PWD/../../programs/OpenMesh-5.1/build/Build/lib/ -lOpenMeshTools

INCLUDEPATH += $$PWD/../../programs/OpenMesh-5.1/build/Build
DEPENDPATH += $$PWD/../../programs/OpenMesh-5.1/build/Build

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/libOpenMeshTools.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/libOpenMeshTools.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/release/OpenMeshTools.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/debug/OpenMeshTools.lib
else:unix: PRE_TARGETDEPS += $$PWD/../../programs/OpenMesh-5.1/build/Build/lib/libOpenMeshTools.a

DISTFILES += \
    extra_code
