
/* This is a place to collect and store compile time configuration settings. */

#define CRKIT_VERSION_MAJOR 1
#define CRKIT_VERSION_MINOR 5
#define CRKIT_VERSION_PATCH 0
#define CRKIT_VERSION_STRING "1.5 ()"

// Remove annoying warning with Visual Studio
#ifdef _MSC_VER
  #pragma warning( disable: 4996 )
#endif

/* #undef Boost_FOUND */
