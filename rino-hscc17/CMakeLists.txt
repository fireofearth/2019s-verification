cmake_minimum_required(VERSION 3.14)
project(RINO_HSCC17)

set(CMAKE_CXX_STANDARD 11)

link_directories(SYSTEM /usr/lib)
link_directories(SYSTEM ~/lib/filib-3.0.2/lib)


include_directories(.)
include_directories(aaflib-0.1)
include_directories(~/lib/fadbad-2.1)
include_directories(~/lib/filib-3.0.2/include)
include_directories( /usr/include)

add_executable(RINO_HSCC17
        aaflib-0.1/aa.h
        aaflib-0.1/aa_aaf.h
        aaflib-0.1/aa_aafapprox.c
        aaflib-0.1/aa_aafarithm.c
        aaflib-0.1/aa_aafcommon.c
        aaflib-0.1/aa_aafinnerflow.c
        aaflib-0.1/aa_aafspecial.c
        aaflib-0.1/aa_aaftrigo.c
        aaflib-0.1/aa_exceptions.c
        aaflib-0.1/aa_exceptions.h
        aaflib-0.1/aa_interval.c
        aaflib-0.1/aa_interval.h
        aaflib-0.1/aa_rounding.c
        aaflib-0.1/aa_rounding.h
        aaflib-0.1/doxygen.cfg
        aaflib-0.1/doxygen.css
        aaflib-0.1/Makefile
        dde_integr.cpp
        dde_integr.h
        fadbad_aa.h
        filib_interval.h
        inner.cpp
        inner.h
        LICENSE
        main.cpp
        Makefile
        matrix.cpp
        matrix.h
        ode_def.cpp
        ode_def.h
        ode_integr.cpp
        ode_integr.h
        README.md
        utils.h)
