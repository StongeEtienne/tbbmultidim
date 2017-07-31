export NTHREADS=10
export TBBDIR=/opt/intel/tbb

#ITKVERSION="4.7.2"

export ITKSettings="-DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF \
-DBUILD_SHARED_LIBS=OFF -DITK_BUILD_DEFAULT_MODULES=OFF \
-DITKGroup_Core=ON -DITKGroup_Numerics=ON -DITKGroup_Filtering=ON \
-DITKGroup_Registration=ON -DModule_ITKIOTransformBase=ON \
-DCMAKE_BUILD_TYPE=Release"
