export SINK_LD_LIBRARY_PATH=/opt/intel/lib/mic/:/opt/intel/tbb/lib/mic/
#/opt/intel/mkl/lib/mic

DATA_PATH=/data

# the data need to be copied on the KNC ***
micnativeloadex build/test_knc/TestFilters -a "$DATA_PATH/nrrd/diffusion.nhdr $DATA_PATH/out/testFilter.nrrd $DATA_PATH/out/testFilter2.nrrd -i 0 -c 8"
micnativeloadex build/test_knc/TestFilters -a "$DATA_PATH/nrrd/diffusion.nhdr $DATA_PATH/out/testFilter.nrrd $DATA_PATH/out/testFilter2.nrrd -i 0 -c 32"
micnativeloadex build/test_knc/TestFilters -a "$DATA_PATH/nrrd/diffusion.nhdr $DATA_PATH/out/testFilter.nrrd $DATA_PATH/out/testFilter2.nrrd -i 0 -c 64"
