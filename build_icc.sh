SCRIPTFOLDER=`dirname $0`
SCRIPTFOLDER=`readlink -f "$SCRIPTFOLDER"`

#OFOLDER=/tmp/itkmic
OFOLDER="${SCRIPTFOLDER}/build/"
SRCDIR="${SCRIPTFOLDER}/src/"

source buildconfig.sh
echo "NTHREADS: $NTHREADS"
echo "TBBDIR: $TBBDIR"

#------------------------------------------------------
# Define subfolders...
#------------------------------------------------------
ITKSOURCE=$OFOLDER/itk
ITKBUILD=$OFOLDER/build_icc
ITKCODE=$OFOLDER/test_icc

PREVDIR=`pwd`
mkdir -p $OFOLDER

#------------------------------------------------------
# First get the source from GIT
#------------------------------------------------------
echo "- Git repository..."
echo "$ITKSOURCE"
if [[ ! -d "$ITKSOURCE" ]]; then
  git clone git://itk.org/ITK.git $ITKSOURCE

  if [[ ! -z "$ITKVERSION" ]]; then
    cd "$ITKSOURCE"
    git checkout "tags/v${ITKVERSION}"
    cd "$PREVDIR"
  fi

else
  cd "$ITKSOURCE"
  svn update
  cd "$PREVDIR"
fi

#------------------------------------------------------
# Configure ITK 
#------------------------------------------------------
lxild=`which xild`
lxiar=`which xiar`
echo
echo "- Configure ..."
mkdir -p $ITKBUILD
echo "$ITKBUILD"
cd $ITKBUILD
CC=icc CXX=icpc CFLAGS=" -axSSE2,SSE4.2,AVX,MIC-AVX512" CXXFLAGS=" -axSSE2,SSE4.2,AVX,MIC-AVX512" \
cmake -DCMAKE_LINKER=${lxild} -DCMAKE_AR=${lxiar} \
-DCFLAGS=" -axSSE2,SSE4.2,AVX,MIC-AVX512" -DCXXFLAGS=" -axSSE2,SSE4.2,AVX,MIC-AVX512" \
$ITKSettings ${ITKSOURCE}/

#------------------------------------------------------
# Build ITK 
#------------------------------------------------------
echo
echo "- Build ..."
make -j $NTHREADS

#------------------------------------------------------
# Configure the sources
#------------------------------------------------------
mkdir -p $ITKCODE
cd "$ITKCODE"
CC=icc CXX=icpc ITK_DIR=$ITKBUILD  CFLAGS=" -lrt -axSSE2,SSE4.2,AVX,MIC-AVX512" CXXFLAGS=" -lrt -std=c++11 -axSSE2,SSE4.2,AVX,MIC-AVX512" \
cmake -DCMAKE_LINKER=${lxild} -DCMAKE_AR=${lxiar} -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_FLAGS=" -lrt -axSSE2,SSE4.2,AVX,MIC-AVX512" -DCMAKE_CXX_FLAGS=" -lrt -std=c++11 -axSSE2,SSE4.2,AVX,MIC-AVX512" \
-DTBB_ROOT=$TBBDIR \
-DITK_DIR=$ITKBUILD ${SRCDIR}/

make -j ${NTHREADS}

