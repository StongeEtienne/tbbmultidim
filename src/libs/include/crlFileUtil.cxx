/**********************************************************************************************//**
 * \file	crlCommon\crlFileUtil.cxx
 *
 * \brief	Useful function to manipulate files 
*************************************************************************************************/

#include "crlFileUtil.h"
#include "crlFileName.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#ifdef WIN32
  #include <windows.h>
  #include <direct.h>
#else
  #include <fnmatch.h>
  #include <dirent.h>
  #include <sys/dir.h>
  #include <unistd.h>
#endif


//-------------------------------------------------
// Create a directory
//-------------------------------------------------

namespace crl {
#ifdef WIN32
#define MKDIR_NONREC(d) mkdir(d)
#else
#define MKDIR_NONREC(d) mkdir(d, S_IRWXU) //0777)
#endif

/**********************************************************************************************//**
 * \fn	void mkdir_rec(const std::string& dirName )
 *
 * \brief	Creates a directory, recursively creating all the subdirectories if needed.
 *
 * \author	Benoit Scherrer
 * \date	September 2015
 *
 * \param	dirName	Pathname of the directory.
 **************************************************************************************************/
void mkdir_rec(const std::string& dirName )
{
	mkdir_rec(dirName.c_str());
}

/**********************************************************************************************//**
 * \fn	void mkdir_rec(const char *dir)
 *
 * \brief	Creates a directory, recursively creating all the subdirectories if needed.
 *
 * \author	Benoit Scherrer
 * \date	September 2015
 *
 * \param	dir	The dir.
 **************************************************************************************************/
void mkdir_rec(const char *dir) 
{
	char tmp[1024];
	char *p = NULL;
	size_t len;

#ifdef WIN32
	_snprintf(tmp, sizeof(tmp),"%s",dir);
#else
	snprintf(tmp, sizeof(tmp),"%s",dir);
#endif
	len = strlen(tmp);
	if(tmp[len - 1] == '/' || tmp[len - 1] == '\\' )
		tmp[len - 1] = 0;
	for(p = tmp + 1; *p; p++)
		if( (*p == '/') || (*p == '\\') ) {
			*p = 0;
			MKDIR_NONREC(tmp);
			*p = '/';
		}

		MKDIR_NONREC(tmp);
}



/**********************************************************************************************//**
 * \fn	bool file_exists(const string& filename)
 *
 * \brief	Check if a file exists. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
 *
 * \param	filename	Filename of the file. 
 *
 * \return	true if the file exists, false if not. 
*************************************************************************************************/
bool file_exists(const string& filename)
{
	struct stat f__stat;
	return (stat(filename.data(),&f__stat) != -1 );
}

bool is_folder(const std::string& dirname )
{
#ifdef WIN32
	DWORD attributs = GetFileAttributes(dirname.c_str());
   	if ( attributs==INVALID_FILE_ATTRIBUTES ) return false;		// does not exists
	return ( (attributs & FILE_ATTRIBUTE_DIRECTORY) == FILE_ATTRIBUTE_DIRECTORY );
#else
	struct stat statbuf;
	if (stat(dirname.data(),&statbuf) == -1 )
		return false;

 	if ( S_ISDIR(statbuf.st_mode) )
		return true;
	else
		return false;

#endif
}

/**********************************************************************************************//**
 * \fn	std::vector<std::string> filelist(const std::string& path, const std::string& pattern )
 *
 * \brief	Returns the list of files in a given directory, for a given pattern.
 *
 * \author	Benoit Scherrer
 * \date	September 2015
 *
 * \param	path   	Full pathname of the file.
 * \param	pattern	The pattern.
 *
 * \return	.
 **************************************************************************************************/
std::vector<std::string> filelist(const std::string& path, const std::string& pattern )
{
	std::vector<std::string> files;

#ifdef WIN32
	//---------------------------------------------
	// Find the first file
	//---------------------------------------------
	WIN32_FIND_DATA ffd;
	std::string str = path+pattern;
	HANDLE hFind = FindFirstFileA(str.c_str(), &ffd);
	if (INVALID_HANDLE_VALUE == hFind) 
		return files;

	//---------------------------------------------
	// Find the next files
	//---------------------------------------------
	do
	{
		if ( !(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			files.push_back(std::string(ffd.cFileName));
		
	}
	while (FindNextFileA(hFind, &ffd) != 0);

	//---------------------------------------------
	// Check if error
	//---------------------------------------------
	DWORD dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES) 
	{
		//Error
	}

	FindClose(hFind);
#else
	

	//std::cout<<"--- SEARCH "<< pattern.c_str() << " in "<< path.c_str()<<std::endl;
	
	//---------------------------------------------
	// Search
	//---------------------------------------------
	DIR *dp;
	struct dirent *ep;
	dp = opendir (path.c_str());
	if (dp == NULL) return files;

	//---------------------------------------------
	// Process all entries
	//---------------------------------------------
	while (ep = readdir (dp))
	{
		//if ( (ep->d_type == DT_DIR) && ( fnmatch(pattern.c_str(), ep->d_name, FNM_FILE_NAME)==0) )
		// DT_DIR doesn t always exists!  and doesn't work on non-ext*	filesystems !
		//  http://comments.gmane.org/gmane.comp.boot-loaders.grub.devel/11855
		if ( !is_folder(path + ep->d_name) && ( fnmatch(pattern.c_str(), ep->d_name, FNM_FILE_NAME)==0) ) 
		{
			//std::cout<<ep->d_name<<std::endl;
			files.push_back(path + ep->d_name);
		}
	}

	//---------------------------------------------
	// Close
	//---------------------------------------------
	closedir (dp);

#endif

	return files;
}



/**********************************************************************************************//**
 * \fn	std::vector<std::string> filelist(const std::string& path, const std::string& pattern )
 *
 * \brief	Returns the list of directories in a given directory, for a given pattern.
 *
 * \author	Benoit Scherrer
 * \date	September 2015
 *
 * \param	path   	Full pathname of the file.
 * \param	pattern	The pattern.
 *
 * \return	.
 **************************************************************************************************/
std::vector<std::string> dirlist(const std::string& path, const std::string& pattern )
{
	std::vector<std::string> files;

#ifdef WIN32
	//---------------------------------------------
	// Find the first file
	//---------------------------------------------
	WIN32_FIND_DATA ffd;
	std::string str = path+pattern;
	HANDLE hFind = FindFirstFileA(str.c_str(), &ffd);
	if (INVALID_HANDLE_VALUE == hFind) 
		return files;

	//---------------------------------------------
	// Find the next files
	//---------------------------------------------
	do
	{
		std::string tmp = std::string(ffd.cFileName);

		if ( (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && (tmp!=".") && (tmp!="..") )
			files.push_back(path + tmp);
		
	}
	while (FindNextFileA(hFind, &ffd) != 0);

	//---------------------------------------------
	// Check if error
	//---------------------------------------------
	DWORD dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES) 
	{
		//Error
	}

	FindClose(hFind);
#else
	

	//std::cout<<"--- SEARCH "<< pattern.c_str() << " in "<< path.c_str()<<std::endl;
	
	//---------------------------------------------
	// Search
	//---------------------------------------------
	DIR *dp;
	struct dirent *ep;
	dp = opendir (path.c_str());
	if (dp == NULL) return files;

	//---------------------------------------------
	// Process all entries
	//---------------------------------------------
	while (ep = readdir (dp))
	{
		//if ( (ep->d_type == DT_DIR) && ( fnmatch(pattern.c_str(), ep->d_name, FNM_FILE_NAME)==0) )
		// DT_DIR doesn t always exists!  and doesn't work on non-ext*	filesystems !
		//  http://comments.gmane.org/gmane.comp.boot-loaders.grub.devel/11855
		if ( (is_folder(path + ep->d_name))  && ( fnmatch(pattern.c_str(), ep->d_name, FNM_FILE_NAME)==0) && (ep->d_name!=".") && (ep->d_name!="..")  ) 
		{
			files.push_back(path + ep->d_name);
		}
	}

	//---------------------------------------------
	// Close
	//---------------------------------------------
	closedir (dp);

#endif

	return files;
}


}


