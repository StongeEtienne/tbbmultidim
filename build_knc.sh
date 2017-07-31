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
ITKBUILD=$OFOLDER/build_knc
ITKCODE=$OFOLDER/test_knc

TOOLCHAIN=$SCRIPTFOLDER/ToolChain.cmake
TRYRUN=$SCRIPTFOLDER/TryRunResults.cmake

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
# Get Intel tools
#------------------------------------------------------
lxild=`which xild`
lxiar=`which xiar`

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

echo "test name print :"
echo "$TOOLCHAIN"
echo "$TRYRUN"
echo "$lxild"
echo "$lxiar"

#------------------------------------------------------
# Configure ITK for knc
#------------------------------------------------------
echo
echo "- Configure ..."
mkdir -p $ITKBUILD
cd $ITKBUILD
CC=icc CXX=icpc \
cmake -DCMAKE_LINKER=${lxild} -DCMAKE_AR=${lxiar} \
$ITKSettings -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN} ${ITKSOURCE}/
#-DCMAKE_CXX_LINK_EXECUTABLE="<CMAKE_AR> -crs  <TARGET> <LINK_FLAGS> <OBJECTS>" \
#-DCMAKE_C_LINK_EXECUTABLE="<CMAKE_AR> -crs  <TARGET> <LINK_FLAGS> <OBJECTS>" \

#------------------------------------------------------
# Configuration with the TryRunResults (ran on MICnativeloadex)
#------------------------------------------------------
echo
echo " -- TryRunResult MicNativeloadex "
echo " copy the file 'TryRunResult.cmake' from $ITKBUILD"
echo " to $OFOLDER/TryRunResult.cmake "
echo " and run on MIC with MicNativeloadex and edit"

cmake -C ${TRYRUN} .

#------------------------------------------------------
# Patch the configuration
#------------------------------------------------------
echo
echo "- Patch SSE2 support ..."
sed -i 's/SSE2_32:INTERNAL=1/SSE2_32:INTERNAL=0/g' $ITKBUILD/CMakeCache.txt
sed -i 's/SSE2_64:INTERNAL=1/SSE2_64:INTERNAL=0/g' $ITKBUILD/CMakeCache.txt

#------------------------------------------------------
# Build ITK for MIC
#------------------------------------------------------
echo
echo "- Build ..."
make -j $NTHREADS


#------------------------------------------------------
# Configure the sources
#------------------------------------------------------
mkdir -p $ITKCODE
cd "$ITKCODE"
CC=icc CXX=icpc ITK_DIR=$ITKBUILD  CFLAGS=" -mmic -tbb " CXXFLAGS=" -mmic -tbb " \
cmake -DCMAKE_LINKER=${lxild} -DCMAKE_AR=${lxiar} -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_FLAGS=" -mmic -tbb " -DCMAKE_CXX_FLAGS=" -mmic -tbb "  \
-DTBB_ROOT=$TBBDIR -DTBB_ARCH_PLATFORM="mic" -DNO_BOOST=1 \
-DQNANHIBIT=1 -HAVE_QNANHIBIT_VALUE=1 -DTEEM_QNANHIBIT=1 \
-DITK_DIR=$ITKBUILD ${SRCDIR}/

make -j ${NTHREADS}

