SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_VERSION 1)
SET(CMAKE_SYSTEM_PROCESSOR k1om)

# specify the cross compiler
SET(CMAKE_C_COMPILER icc)
SET(CMAKE_CXX_COMPILER icpc)


# archive and linker (DOESNT WORK IF SET AFTER THE CMAKE CALL)
# NEED TO BE SET IN THE COMMAND LINE
#SET(CMAKE_AR "/opt/intel/composer_xe_2015.3.187/bin/intel64/xiar")
#SET(CMAKE_LINKER "/opt/intel/composer_xe_2015.3.187/bin/intel64/xild")
#SET(CMAKE_C_LINK_EXECUTABLE "<CMAKE_AR> -crs <TARGET> <LINK_FLAGS> <OBJECTS>")
#SET(CMAKE_CXX_LINK_EXECUTABLE "<CMAKE_AR> -crs <TARGET> <LINK_FLAGS> <OBJECTS>")

# cmake flag
#SET(CMAKE_C_FLAGS "-mmic")
#SET(CMAKE_CXX_FLAGS "-mmic")
SET(CMAKE_C_FLAGS "-mmic " CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS "-mmic -std=c++11 " CACHE STRING "" FORCE)



#SET(MPI_C_COMPILER mpiicc)
SET(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)#todo change

# ITK specific
# where is the target environment
#SET(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)#todo change
SET(CMAKE_FIND_ROOT_PATH
  "/opt/mpss/3.5.2/sysroots/x86_64-mpsssdk-linux/"
  "/opt/mpss/3.6.1/sysroots/x86_64-mpsssdk-linux/"
  "/usr/linux-k1om-4.7/"
)
# search for programs in the build host directories
#SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
#SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# where is the target environment 

