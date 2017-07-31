/**********************************************************************************************//**
 * \file	crlLogClass.cxx
 *
 * \brief	The crl::LogClass class. 
*************************************************************************************************/

/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "configuration.h"
#include "crlLogClass.h"
#include "crlFileName.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <time.h>

#include "itkTimeProbe.h"

#ifdef WIN32
  #include <windows.h>
  #pragma warning( disable: 4996 )
#endif

using namespace std;

namespace crl {		


#ifdef WIN32
    #define LOCK_MUTEX WaitForSingleObject(d->lockMutex, INFINITE);
    #define UNLOCK_MUTEX ReleaseMutex(d->lockMutex);
#else
    #define LOCK_MUTEX pthread_mutex_lock(&d->lockMutex);
    #define UNLOCK_MUTEX pthread_mutex_unlock(&d->lockMutex);
#endif


	/**********************************************************************************************//**
	 * \class	LogClass::PrivateData
	 *
	 * \brief	The private data for LogClass. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	class LogClass::PrivateData
	{
	public:
		PrivateData() {}

		std::fstream	ofile;

		clock_t			timerStart, timerEnd ;
		time_t			timerRealStart, timerRealEnd ;

		itk::TimeProbe	timerItk;

#ifdef WIN32
		HANDLE  lockMutex;
#else
		pthread_mutex_t lockMutex;
#endif

	};

	/**********************************************************************************************//**
	 * \fn	LogClass::LogClass()
	 *
	 * \brief	Constructor.
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 **************************************************************************************************/
	LogClass::LogClass()
	{
		d = new PrivateData;

#ifdef WIN32
		d->lockMutex = CreateMutex(NULL, FALSE, NULL);
#else
		pthread_mutex_init(&d->lockMutex, NULL);
#endif
	}

	/**********************************************************************************************//**
	 * \fn	LogClass::~LogClass()
	 *
	 * \brief	Destructor.
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 **************************************************************************************************/
	LogClass::~LogClass()
	{
		if ( d->ofile.is_open() ) d->ofile.close();
		delete d;
	}

	/**********************************************************************************************//**
	 * \fn	void LogClass::setOutputBaseFileName(const std::string& baseFileName)
	 *
	 * \brief	Sets the output base file name. The log file will be the output base file name with
	 * 			the extension '.log'.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \param	baseFileName	Filename of the base file.
	 **************************************************************************************************/
	void LogClass::setOutputBaseFileName(const std::string& baseFileName)
	{
		LOCK_MUTEX;
		crl::FileName fileName(baseFileName);
		std::string fn = fileName.getPath()+fileName.getFileName()+".log";
		d->ofile.open ( fn.c_str(), std::fstream::out);
		if  ( !d->ofile.is_open() ) 
		{
			std::cout<<"- WARNING. Cannot open <" << fn << "> for writting."<<std::endl;
		}
		UNLOCK_MUTEX;
	}

	/**********************************************************************************************//**
	 * \fn	void LogClass::flush()
	 *
	 * \brief	Flushs the output log file.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 **************************************************************************************************/
	void LogClass::flush()
	{
		LOCK_MUTEX;
		d->ofile.flush();
		UNLOCK_MUTEX;
	}

	/**********************************************************************************************//**
	 * \fn	void LogClass::close()
	 *
	 * \brief	Closes the output log file.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 **************************************************************************************************/
	void LogClass::close()
	{
		LOCK_MUTEX;
		d->ofile.close();
		UNLOCK_MUTEX;
	}

	void LogClass::outputDefaultLogHeader(int argc, char **argv, const std::string& toolname, const std::string& copyright, const std::string& toolversion )
	{
		LOCK_MUTEX;
		if ( d->ofile.is_open() )
		{
			d->ofile<<"================================================="<<endl;
			d->ofile<<" " << toolname << endl;
			if ( toolversion!="" ) d->ofile<<" " << toolversion << endl;
			d->ofile << std::endl;
			if ( copyright!="" ) d->ofile<<" " << copyright << endl;
			d->ofile<<" CRKIT VERSION: "<<CRKIT_VERSION_STRING<<endl;
			d->ofile<<"-------------------------------------------------"<<endl;
			d->ofile<<" STARTED ON " <<currentDateTime()<<std::endl;
			d->ofile<<"================================================="<<endl<<endl;
			d->ofile<<"COMMAND LINE"<<std::endl;
			d->ofile<<"------------"<<std::endl;
			d->ofile<<cmdLine(argc,argv)<<std::endl<<std::endl;
			d->ofile.flush();
		}
		UNLOCK_MUTEX;
	}

	/**********************************************************************************************//**
	 * \fn	std::string LogClass::cmdLine(int argc, char **argv)
	 *
	 * \brief	Gets the command line as a std::string
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \param	argc			The argc.
	 * \param [in,out]	argv	If non-null, the argv.
	 *
	 * \return	.
	 **************************************************************************************************/
	std::string LogClass::cmdLine(int argc, char **argv)
	{
		std::string str;

		for ( int i=0; i<argc; i++ )
		{
			if ( i!=0 ) str += " ";
			str += std::string(argv[i]);
		}
		return str;
	}

	/**********************************************************************************************//**
	 * \fn	std::string LogClass::currentDateTime() const
	 *
	 * \brief	Gets the current date and time as a std::string.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \return	.
	 **************************************************************************************************/
	std::string LogClass::currentDateTime() const
	{
		time_t		currentTime;
		struct tm	*timeinfo;

		time ( &currentTime );
		timeinfo = localtime ( &currentTime );
		return std::string( asctime (timeinfo) );
	}

	/**********************************************************************************************//**
	 * \fn	void LogClass::timerStart()
	 *
	 * \brief	Starts the internal timer to measure time.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 **************************************************************************************************/
	void LogClass::timerStart()
	{
		LOCK_MUTEX;
		d->timerStart  = clock();  
		d->timerEnd = 0;

		d->timerRealStart = time(NULL);
		d->timerRealEnd = time(NULL); 

		d->timerItk.Start();
		UNLOCK_MUTEX;
	}

	/**********************************************************************************************//**
	 * \fn	void LogClass::timerStop()
	 *
	 * \brief	Stops the internal timer
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 **************************************************************************************************/
	void LogClass::timerStop()
	{
		LOCK_MUTEX;
		d->timerEnd = clock(); 
		d->timerRealEnd = time(NULL);
		d->timerItk.Stop();
		UNLOCK_MUTEX;
	}

	/**********************************************************************************************//**
	 * \fn	float LogClass::timerRealProcessorLaps() const
	 *
	 * \brief	Gets the timer real laps.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \return	.
	 **************************************************************************************************/
	float LogClass::timerRealProcessorLaps() const
	{
		LOCK_MUTEX;
		float delay= (float)(d->timerEnd - d->timerStart);
		//float delay= (float)difftime(d->timerRealEnd, d->timerRealStart); // no not this one
		delay/=CLOCKS_PER_SEC;// divide by the number of clock ticks per second // (delay/=CLK_TCK) on Borland Windows?
		UNLOCK_MUTEX;

		return(delay);
	}

	/**********************************************************************************************//**
	 * \fn	float LogClass::timerLaps() const
	 *
	 * \brief	Gets the timer laps.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \return	.
	 **************************************************************************************************/
	float LogClass::timerLaps() const
	{
		LOCK_MUTEX;
		float delay = d->timerItk.GetTotal();
		UNLOCK_MUTEX;

		return(delay);
	}

	/**********************************************************************************************//**
	 * \fn	std::string LogClass::timerFormatRealProcessorLaps() const
	 *
	 * \brief	Gets the timer real laps and format it as a string hh::mm::ss.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \return	.
	 **************************************************************************************************/
	std::string	 LogClass::timerFormatRealProcessorLaps() const
	{
		return formatTime(timerRealProcessorLaps());
	}

	/**********************************************************************************************//**
	 * \fn	std::string LogClass::timerFormatLaps() const
	 *
	 * \brief	Gets the timer format laps and format it as a string hh::mm::ss.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \return	.
	 **************************************************************************************************/
	std::string LogClass::timerFormatLaps() const
	{
		return formatTime(timerLaps());
	}

	/**********************************************************************************************//**
	 * \fn	std::string LogClass::formatTime(float seconds) const
	 *
	 * \brief	Format a time (hh::mm::ss) given in seconds.
	 *
	 * \author	Benoit Scherrer
	 * \date	December 2011
	 *
	 * \param	seconds	The seconds.
	 *
	 * \return	The formatted time.
	 **************************************************************************************************/
	std::string LogClass::formatTime(float seconds) const
	{
		int h, mn, s;
		char szBuffer[128];

		h = (int)(seconds/3600.0f);  seconds -= 3600.0f*h;
		mn = (int)(seconds/60.0f);   seconds -= 60.0f*mn;
		s = (int)(seconds+0.5f); 

		if ( h>0 )
			sprintf(szBuffer,"%02dh%02dm%02ds", h, mn, s);
		else if (mn>0 )
			sprintf(szBuffer,"%02dm%02ds", mn, s);
		else
			sprintf(szBuffer,"%02ds", s);

		return std::string(szBuffer);
	}

	/**********************************************************************************************//**
	 * \property	std::ostream& LogClass::operator<<(const std::string& str)
	 *
	 * \brief	<< operator
	 *
	 * \return	
	 **************************************************************************************************/
	std::ostream& LogClass::operator<<(const std::string& str)
	{
		LOCK_MUTEX;
		if ( d->ofile.is_open() )
		{
			d->ofile << str;
		}
		UNLOCK_MUTEX;
		return d->ofile;
		
	}


	std::ostream& LogClass::operator<<(int i)
	{
		LOCK_MUTEX;
		if ( d->ofile.is_open() )
		{
			d->ofile << i;
		}
		UNLOCK_MUTEX;
		return d->ofile;
	}

	std::ostream& LogClass::operator<<(float f)
	{
		LOCK_MUTEX;
		if ( d->ofile.is_open() )
		{
			d->ofile << f;
		}
		UNLOCK_MUTEX;
		return d->ofile;
	}
}
