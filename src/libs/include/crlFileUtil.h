/**********************************************************************************************//**
 * \file	crlCommon\crlFileUtil.h
 *
 * \brief	Declares useful functions to manipulate files 
*************************************************************************************************/

#ifndef h_CRL_FILE_UTIL
#define h_CRL_FILE_UTIL

#include <string>
#include <vector>
using namespace std;

namespace crl {
	extern void mkdir_rec( const std::string& dirName );
	extern void mkdir_rec( const char *dir );

	extern bool file_exists(const string& filename);
        extern bool is_folder(const std::string& dirname );

	extern std::vector<std::string> filelist(const std::string& path, const std::string& pattern );
	extern std::vector<std::string> dirlist(const std::string& path, const std::string& pattern );
}

#endif
