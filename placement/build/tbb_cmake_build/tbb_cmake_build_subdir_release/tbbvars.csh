#!/bin/csh
setenv TBBROOT "/home/zec022@AD.UCSD.EDU/placement/build/oneTBB-2019_U9" #
setenv tbb_bin "/home/zec022@AD.UCSD.EDU/placement/build/tbb_cmake_build/tbb_cmake_build_subdir_release" #
if (! $?CPATH) then #
    setenv CPATH "${TBBROOT}/include" #
else #
    setenv CPATH "${TBBROOT}/include:$CPATH" #
endif #
if (! $?LIBRARY_PATH) then #
    setenv LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LIBRARY_PATH "${tbb_bin}:$LIBRARY_PATH" #
endif #
if (! $?LD_LIBRARY_PATH) then #
    setenv LD_LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LD_LIBRARY_PATH "${tbb_bin}:$LD_LIBRARY_PATH" #
endif #
 #
