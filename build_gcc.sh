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
ITKBUILD=$OFOLDER/build_gcc
ITKCODE=$OFOLDER/test_gcc

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
lxild=`which ld`
lxiar=`which ar`
echo
echo "- Configure ..."
mkdir -p $ITKBUILD
echo "$ITKBUILD"
cd $ITKBUILD
cmake -DCMAKE_LINKER=${lxild} -DCMAKE_AR=${lxiar} $ITKSettings ${ITKSOURCE}/

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
ITK_DIR=$ITKBUILD   CXXFLAGS=" -std=c++11 " \
cmake -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS=" -std=c++11 " \
-DTBB_ROOT=$TBBDIR \
-DITK_DIR=$ITKBUILD ${SRCDIR}/

make -j ${NTHREADS}

