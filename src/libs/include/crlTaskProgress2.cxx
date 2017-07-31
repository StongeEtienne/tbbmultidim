/**********************************************************************************************//**
 * \file	crlTaskProgress2.cxx
 *
 * \brief	The crl::TaskProgress2 class. 
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


#include "crlTaskProgress2.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <cstdio>

#ifdef WIN32
  #include <windows.h>
  #pragma warning( disable: 4996 )
#endif

using namespace std;

namespace crl {

/*--------------------------------------------------------------
 Some VT100 macros (don t work for win32)
--------------------------------------------------------------*/
#define vt100_reset "\33c"
#define vt100_clear_screen "\33[2J"
#define vt100_cursor_up(    count ) "\33[" << count << "A"
#define vt100_cursor_down(  count ) "\33[" << count << "B"
#define vt100_cursor_right( count ) "\33[" << count << "C"
#define vt100_cursor_left(  count ) "\33[" << count << "D"

#define vt100_cursor_savepos "\33[s"
#define vt100_cursor_restorepos "\33[u"


/**********************************************************************************************//**
 * \class	TaskProgress2::PrivateData
 *
 * \brief	The private data for FileName. 
 *
 * \author	Benoit Scherrer
 * \date	January 2010
*************************************************************************************************/
class TaskProgress2::PrivateData
	{
	public:
		PrivateData( unsigned long int counterMaxVal, bool showRemainingTime ):
			m_CounterMaxVal(counterMaxVal),
			m_CurrentCounter(0),
			m_CounterBeforeUpdate((long int)(m_CounterMaxVal/100 + 0.5)),
			m_Silent(false),
			m_ShowRemainingTime(showRemainingTime),
			m_LinePrefix("")
			{}

		itk::TimeProbe				m_TimeProbe;
		itk::SimpleFastMutexLock	m_ProgressMutex;


		unsigned long int	m_CounterMaxVal;
		unsigned long int	m_CurrentCounter;
		long int			m_CounterBeforeUpdate;
		std::vector<int>	m_TimeFilter;

		bool				m_Silent;
		bool				m_ShowRemainingTime;
		std::string			m_LinePrefix;
		std::string			m_LineSuffix;
	};

	/**********************************************************************************************//**
	 * \fn	TaskProgress2::TaskProgress2( unsigned long int counterMaxVal, bool showRemainingTime )
	 *
	 * \brief	Constructor.
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	counterMaxVal	 	The maximum value of the counter (corresponding to 100%).
	 * \param	showRemainingTime	If true, will show the estimated remaining time.
	 **************************************************************************************************/
	TaskProgress2::TaskProgress2( unsigned long int counterMaxVal, bool showRemainingTime )
	{
		d = new TaskProgress2::PrivateData( counterMaxVal, showRemainingTime );
	}

	/**********************************************************************************************//**
	 * \fn	TaskProgress2::~TaskProgress2()
	 *
	 * \brief	Destructor. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	TaskProgress2::~TaskProgress2()
	{
		delete d;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::InitializeProgress(unsigned long int counterMaxVal, bool startTimer )
	 *
	 * \brief	Initializes the progress. Provides the maximum counter value to reach (=100%), and
	 * 			start the progress if needed.
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2012
	 *
	 * \param	counterMaxVal	If not null, set the counter maximum value. Else will keep the
	 * 							previously set counter maximum value.
	 * \param	startTimer   	true to start timer.
	 **************************************************************************************************/
	void TaskProgress2::InitializeProgress(unsigned long int counterMaxVal, bool startTimer )
	{
		d->m_TimeProbe = itk::TimeProbe();
		d->m_CurrentCounter = 0;
		d->m_CounterBeforeUpdate = 0;
		d->m_TimeFilter.clear();
		if ( counterMaxVal!=0 ) d->m_CounterMaxVal = counterMaxVal;

		if ( startTimer ) StartProgress();
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::StdOutLock()
	 *
	 * \brief	Mutex mechanism. Lock any output on std out (be careful. Must be unlocked, else
	 * 			risk of deadlock!)
	 *
	 * \author	Benoit Scherrer
	 * \date	November 2014
	 **************************************************************************************************/
	void TaskProgress2::StdOutLock()
	{
		d->m_ProgressMutex.Lock();
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::StdOutUnlock()
	 *
	 * \brief	Mutex mechanism. Unlock the output on std out
	 *
	 * \author	Benoit Scherrer
	 * \date	November 2014
	 **************************************************************************************************/
	void TaskProgress2::StdOutUnlock()
	{
		d->m_ProgressMutex.Unlock();
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::StartProgress()
	 *
	 * \brief	Starts the progress.
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2012
	 **************************************************************************************************/
	void TaskProgress2::StartProgress()
	{
		d->m_TimeProbe.Start();

		d->m_ProgressMutex.Lock();
		ShowProgress();
		d->m_ProgressMutex.Unlock();
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::Update()
	 *
	 * \brief	Increment the internal counter and print out the result if needed
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	*************************************************************************************************/
	void TaskProgress2::IncrementProgress()
	{
		d->m_ProgressMutex.Lock();

		d->m_CounterBeforeUpdate--;
		d->m_CurrentCounter++;

		if ( d->m_CounterBeforeUpdate<=0 || d->m_CurrentCounter==d->m_CounterMaxVal ) 
		{
			d->m_CounterBeforeUpdate = (long int)(d->m_CounterMaxVal/100 + 0.5);
			ShowProgress();
		}
                d->m_ProgressMutex.Unlock();
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::ShowProgress()
	 *
	 * \brief	Internal function. Shows the current progress. WARNING: d->m_ProgressMutex.Lock();
	 * 			must have been called before calling this function!
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2012
	 **************************************************************************************************/
	void TaskProgress2::ShowProgress()
	{
		if ( d->m_Silent ) return;

#ifdef WIN32
		//-------------------------------------
		// If we are under win32, we can use 
		//functions to get/set the cursor pos
		//-------------------------------------
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		COORD pos;
		if ( GetConsoleScreenBufferInfo( GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
			pos.Y = csbi.dwCursorPosition.Y;
		else
			pos.Y = 0;
		pos.X=0;
		SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);
#else
		//-------------------------------------
		// Linux/OSX: save current cursor position
		//-------------------------------------
		//std::cout<<vt100_cursor_left(40); 
		//std::cout<<vt100_cursor_up(1);
		std::cout << vt100_cursor_savepos ;
#endif

		//-------------------------------------
		// Get the percentage of the task completed
		//-------------------------------------
		float percentage ;
		if ( d->m_CounterMaxVal==0 ) percentage=0.0f;
		else percentage = (float) ( d->m_CurrentCounter * 100.0 / ((double) d->m_CounterMaxVal ));

		//-------------------------------------
		// Convert to a string
		//-------------------------------------
		char szBuffer[128];
		sprintf(szBuffer,"%03d", (int)percentage);

		//-------------------------------------
		// Check that we do not over count
		//-------------------------------------
		if ( d->m_CurrentCounter <= d->m_CounterMaxVal ) 
		{
			std::cout<< d->m_LinePrefix.c_str() << szBuffer <<"%  "; 

			//-------------------------------------
			// If the task is finished (100%), show the
			// total time
			//-------------------------------------
			if ( d->m_CurrentCounter==d->m_CounterMaxVal ) 
			{
				d->m_TimeProbe.Stop();
				int t = (int)d->m_TimeProbe.GetTotal();
				std::string strTime = ConvertTimeToString(t);
				std::cout<< "("<<strTime <<") " <<d->m_LineSuffix.c_str() << "                          ";
				#ifdef WIN32
					std::cout<<std::endl;
				#endif
			}
			//-------------------------------------	
			// Else show the estimated left time 
			//-------------------------------------
			else 
			{
				d->m_TimeProbe.Stop();
				itk::TimeProbe::TimeStampType time = d->m_TimeProbe.GetTotal();
				d->m_TimeProbe.Start();

				//std::cout<< "c="<<d->m_CurrentCounter<< " max=" << d->m_CounterMaxVal<< "t=" << time;

				if ( percentage!=0.0f && d->m_ShowRemainingTime ) 
				{
					int predictedTime = (int)(100.0*time/percentage);

					if ( predictedTime-time != 0 )
					{
						//-----------------------------------------
						// Compute a filtered estimated time left
						//-----------------------------------------
						double filteredPredictedTime=0;

						if ( d->m_TimeFilter.size()>5 ) d->m_TimeFilter.pop_back();
						d->m_TimeFilter.insert(d->m_TimeFilter.begin(), predictedTime);

						for ( unsigned int i=0; i<d->m_TimeFilter.size(); i++ )
							filteredPredictedTime+=d->m_TimeFilter[i];
						filteredPredictedTime=filteredPredictedTime/d->m_TimeFilter.size();

						//std::cout<<" "<<predictedTime<<" - " << filteredPredictedTime << "  ";

						std::string strTime = ConvertTimeToString(filteredPredictedTime-time);
						std::cout<< "(Estimated left: "<<strTime.c_str()<<") " << d->m_LineSuffix.c_str() << "               ";
					}
				}

				if ( !d->m_ShowRemainingTime ) std::cout<<d->m_LineSuffix.c_str();
			}

			std::flush(std::cout);
		}
#ifndef WIN32
		std::cout << vt100_cursor_restorepos ;
#endif

	}


	/**********************************************************************************************//**
	 * \fn	std::string GetTotalTime() const
	 *
	 * \brief	Gets the total time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	April 2010
	 *
	 * \return	The total time. 
	*************************************************************************************************/
	std::string TaskProgress2::ConvertTimeToString(itk::RealTimeClock::TimeStampType t) 
	{
		char szBuffer[512];
	
		int ss = ((int)t)%60; 
		int mm = (int)((t-ss)/60)%60;
		int hh = (int)(t / 3600);
		sprintf(szBuffer, "%02dh %02dm %02ds", hh, mm, ss );

		return string(szBuffer);
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::SetLinePrefix(const std::string& prefix)
	 *
	 * \brief	Sets a line prefix (used for cout). 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	prefix	The prefix. 
	*************************************************************************************************/
	void TaskProgress2::SetLinePrefix(const std::string& prefix)
	{
		d->m_LinePrefix = prefix;
	}

	/**********************************************************************************************//**
	 * \fn	std::string TaskProgress2::GetLinePrefix()
	 *
	 * \brief	Gets the line prefix. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \return	The line prefix. 
	*************************************************************************************************/
	std::string TaskProgress2::GetLinePrefix() const
	{
		return d->m_LinePrefix;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::SetLinePrefix(const std::string& prefix)
	 *
	 * \brief	Sets a line prefix (used for cout). 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	prefix	The prefix. 
	*************************************************************************************************/
	void TaskProgress2::SetLineSuffix(const std::string& suffix)
	{
		d->m_LineSuffix = suffix;
	}

	/**********************************************************************************************//**
	 * \fn	std::string TaskProgress2::GetLinePrefix()
	 *
	 * \brief	Gets the line prefix. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \return	The line prefix. 
	*************************************************************************************************/
	std::string TaskProgress2::GetLineSuffix() const
	{
		return d->m_LineSuffix;
	}


	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::SetShowEstimatedRemainingTime(bool b)
	 *
	 * \brief	Sets if it should show the estimated remaining time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \param	b	true to. 
	*************************************************************************************************/
	void TaskProgress2::SetShowEstimatedRemainingTime(bool b)
	{
		d->m_ShowRemainingTime = b;
	}

	/**********************************************************************************************//**
	 * \fn	bool TaskProgress2::GetShowEstimatedRemainingTime()
	 *
	 * \brief	Gets if it shows the estimated remaining time. 
	 *
	 * \author	Benoit Scherrer
	 * \date	January 2010
	 *
	 * \return	true if it succeeds, false if it fails. 
	*************************************************************************************************/
	bool TaskProgress2::GetShowEstimatedRemainingTime() const
	{
		return d->m_ShowRemainingTime;
	}

	/**********************************************************************************************//**
	 * \fn	void TaskProgress2::SetSilent(bool b)
	 *
	 * \brief	Sets if silent. Can be usefull for some algorithms to have a silent mode
	 *
	 * \author	Benoit Scherrer
	 * \date	November 2013
	 *
	 * \param	b	true to b.
	 **************************************************************************************************/
	void TaskProgress2::SetSilent(bool b) 
	{
		d->m_Silent = b;
	}

	/**********************************************************************************************//**
	 * \fn	bool TaskProgress2::GetSilent() const
	 *
	 * \brief	Gets if silent.
	 *
	 * \author	Benoit Scherrer
	 * \date	November 2013
	 *
	 * \return	true if it succeeds, false if it fails.
	 **************************************************************************************************/
	bool  TaskProgress2::GetSilent() const
	{
		return d->m_Silent;
	}






	//ItkFilterProgressObserver::ItkFilterProgressObserver():
	//m_TaskProgress2(0,1,100,true, false) 
	//{ }

	//void ItkFilterProgressObserver::Execute(itk::Object *caller, const itk::EventObject & event)
	//{
	//	Execute( (const itk::Object *)caller, event);
	//}

	//void ItkFilterProgressObserver::Execute(const itk::Object * object, const itk::EventObject & event)
	//{
	//	ProcessPointer filter = dynamic_cast< ProcessPointer >( object );
	//	if( typeid( event ) == typeid( itk::ProgressEvent ) )
	//		m_TaskProgress2.SetCurrentCounter((unsigned int)(filter->GetProgress()*100));

	//	if( typeid( event ) == typeid( itk::StartEvent ) )
	//	{
	//		TaskProgress2::InitTaskProgress2();
	//		m_TaskProgress2.SetCurrentCounter(0);
	//	}
	//	//			if( typeid( event ) == typeid( itk::EndEvent ) )
	//}

	//void ItkFilterProgressObserver::SetProgressLinePrefix(const std::string& str)
	//{
	//	m_TaskProgress2.SetLinePrefix(str);
	//}
}
