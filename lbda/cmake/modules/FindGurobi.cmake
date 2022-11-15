include(FindPackageHandleStandardArgs)

set(GUROBI_HOME "$ENV{GUROBI_HOME}")

find_library(
        GUROBI_C_LIBRARY
        NAMES libgurobi90.so gurobi90.lib libgurobi100.so gurobi100.lib #added version 10 names        
        PATHS "${GUROBI_HOME}/lib")

find_library(GUROBI_CXX_LIBRARY
        NAMES libgurobi_c++.a gurobi_c++*.lib gurobi_c++md2017 #should be ok, added one anyway
        PATHS "${GUROBI_HOME}/lib")

set(GUROBI_INCLUDES ${GUROBI_HOME}/include)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(GUROBI DEFAULT_MSG
        GUROBI_INCLUDES
        GUROBI_CXX_LIBRARY
        GUROBI_C_LIBRARY)

if (GUROBI_FOUND)
    set(GUROBI_INCLUDE_DIRS ${GUROBI_HOME}/include)
    set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_C_LIBRARY}")

    mark_as_advanced(
            GUROBI_CXX_LIBRARY
            GUROBI_C_LIBRARY
            GUROBI_HOME)
endif ()
