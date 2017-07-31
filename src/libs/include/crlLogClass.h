/**********************************************************************************************//**
 * \file	crlLogClass.h
 *
 * \brief	Declares the crl::LogClass class. 
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


#ifndef h_CRL_LOG_CLASS
#define h_CRL_LOG_CLASS

#ifdef WIN32
#else
#include <pthread.h>
#endif

#include <string>
#include <itkTimeProbe.h>
#include <itkCommand.h>
#include <itkProcessObject.h>


namespace crl {

	// need to include "tclap/CmdLine.h" in your main to use that  (#include "tclap/CmdLine.h")
    #define TCLAP_LOG_SWITCH(cmd) TCLAP::SwitchArg argLogSwitch("","log", "Output log information in [output basefile].log.", cmd, false);


	/**********************************************************************************************//**
	 * \class	LogClass
	 *
	 * \brief	
	 *
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	class LogClass
	{
	public:
		LogClass();
		~LogClass();

		void			setOutputBaseFileName(const std::string& baseFileName);
		void			flush();
		void			close();

		void			outputDefaultLogHeader(int argc, char **argv, const std::string& toolname, const std::string& copyright="", const std::string& toolversion="" );

		std::string		cmdLine(int argc, char **argv);

		std::string		currentDateTime() const;

		void			timerStart();
		void			timerStop();
		float			timerRealProcessorLaps() const;
		float			timerLaps() const;
		std::string		timerFormatRealProcessorLaps() const;
		std::string		timerFormatLaps() const;

		std::string		formatTime(float seconds) const;

		std::ostream& operator<<(const std::string& str);
		std::ostream& operator<<(int i);
		std::ostream& operator<<(float f);

	protected:

	private:
		class PrivateData;
		PrivateData *d;

	};

	

}

#endif
