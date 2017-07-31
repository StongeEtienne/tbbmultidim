===============================================================

  crlMFMEstimate Bundle for Intel

  2017/06/26

  (c) Benoit Scherrer, Etienne St-Onge, 2017
      Computational Radiology Laboratory (CRL)
      Boston Children's Hospital

===============================================================
___________________________________

# SETUP :

Type of build: (look below for more details)
----------------------------------------
    I - Standart compilation
       a) GCC compiler
            ./build_gcc.sh
            ./test_gcc.sh
      ----------------------------------
       b) Intel compiler
            ./build_icc.sh
            ./test_icc.sh
   -------------------------------------
    II - MIC (Intel compiler)
       a) Knight Landing (KNL)
            ./build_icc.sh
            ./test_icc.sh
      ----------------------------------
       b) Knight Corner (KNC)
            #copy data to KNC (test_knc.sh DATA_PATH)
            ./build_knc.sh
            ./test_knc.sh

----------------------------------------

Variable that may need to be changed in the build script (buildconfig.sh): 

    # TBB path:
    TBBDIR=/opt/intel/tbb

    # Number of thread to compile:
    NTHREADS=8

    # Debug:
    change DCMAKE_BUILD_TYPE=Release
    to     DCMAKE_BUILD_TYPE=ReleaseWithDebugInfo
    or     DCMAKE_BUILD_TYPE=Debug


Parameters for benchmark (test_*.sh):
    # -c NB_CORES
    # -i method for multithreading 
      (0: test all, 1:itk, 2:itk multidim, 3:tbb, 4:tbb multidim)

_______________________________________________________________
_______________________________________________________________
# Compilation debug with icc
-----------------
    # Might be required with intel compiler
        source /opt/intel/bin/compilervars.sh intel64
        source /opt/intel/tbb/bin/tbbvars.sh  intel64

