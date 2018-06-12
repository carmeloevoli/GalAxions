cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(ymw16 STATIC
        dora.c
        fermibubble.c
        galcen.c
        gum.c
        localbubble.c
        lmc.c
        nps.c
        smc.c
        spiral.c
        thick.c
        thin.c
        ymw16.c
        ymw16par.c
        ne.c
        ne_crd.c
        )

SET_TARGET_PROPERTIES(ymw16 PROPERTIES COMPILE_FLAGS "-O2")

