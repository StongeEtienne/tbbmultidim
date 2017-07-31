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
#ifndef __crlMultiThreadedWorkPileAbstract_txx
#define __crlMultiThreadedWorkPileAbstract_txx

#include "crlMultiThreadedWorkPileAbstract.h"
#include <math.h>
#include <stdio.h>

/*--------------------------------------------------------------
 Some VT100 macros (don t work for win32)
--------------------------------------------------------------*/
#define vt100_clear_screen "\33[2J"
#define vt100_cursor_up(    count ) "\33[" << count << "A"
#define vt100_cursor_down(  count ) "\33[" << count << "B"
#define vt100_cursor_right( count ) "\33[" << count << "C"
#define vt100_cursor_left(  count ) "\33[" << count << "D"

#define vt100_cursor_savepos "\33[s"
#define vt100_cursor_restorepos "\33[u"

namespace crl {

/**********************************************************************************************//**
 * \fn	template< class TJobData > crlMultiThreadedWorkPile<TJobData>::crlMultiThreadedWorkPile()
 *
 * \brief	Default constructor.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \tparam	TJobData	Type of the work pile job.
 **************************************************************************************************/
template< class TJobData >
crlMultiThreadedWorkPileAbstract<TJobData>::crlMultiThreadedWorkPileAbstract()
{
	this->m_MultiThreader       = itk::MultiThreader::New();
	this->m_ProgressTimeProbe	= itk::TimeProbe();
	this->m_ShowProgress		= true;
	this->m_NbThreads			= 1;
	this->m_UseAutomaticProgressCounter = true;

}

/**********************************************************************************************//**
 * \fn	template< class TJobData > crlMultiThreadedWorkPile<TJobData>::::~WindowedMultiThreadedImageFilter()
 *
 * \brief	Destructor.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \tparam	TJobData	Type of the work pile job.
 **************************************************************************************************/
template< class TJobData >
crlMultiThreadedWorkPileAbstract<TJobData>::~crlMultiThreadedWorkPileAbstract()
{

}

/**********************************************************************************************//**
 * \fn	template< class TInputImage, class TOutputImage > itk::MultiThreader* WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::GetMultiThreader() const
 *
 * \brief	Gets the multi threader.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \return	The multi threader.
 **************************************************************************************************/
template< class TJobData >
itk::MultiThreader* crlMultiThreadedWorkPileAbstract<TJobData>::GetMultiThreader() const
{
	return this->m_MultiThreader;
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPileAbstract<TJobData>::LockJobQueueAccess()
 *
 * \brief	Locks the job queue access (mutex).
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::LockJobQueueAccess()
{
	this->m_JobQueueMutex.Lock();
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPileAbstract<TJobData>::UnlockJobQueueAccess()
 *
 * \brief	Unlocks the job queue access (mutex).
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::UnlockJobQueueAccess()
{
	this->m_JobQueueMutex.Unlock();
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPileAbstract<TJobData>::ExecuteQueue()
 *
 * \brief	Executes all the jobs in the queue.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::ExecuteQueue()
{
	//-------------------------------------
	// Create the barrier
	//-------------------------------------
	this->m_StartBarrierSync = itk::Barrier::New();
	this->m_StartBarrierSync->Initialize(this->m_NbThreads);
	this->m_FinishBarrierSync = itk::Barrier::New();
	this->m_FinishBarrierSync->Initialize(this->m_NbThreads);

	//-------------------------------------
	// Prepare the multi-threader
	//-------------------------------------
	this->m_MultiThreader->SetNumberOfThreads( this->m_NbThreads );
	this->m_MultiThreader->SetSingleMethod( this->WorkPileThreaderCallback, (void *)this );

	//-------------------------------------
	// Start progress time probe
	//------------------------------------- 
	if ( this->GetShowProgress() ) 
		this->StartProgressTimeProbe();

	//-------------------------------------
	// Run!
	//-------------------------------------
	this->m_MultiThreader->SingleMethodExecute();
}

/**********************************************************************************************//**
 * \fn	template < class TInputImage, class TOutputImage > ITK_THREAD_RETURN_TYPE WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::WindowedThreaderCallback( void *arg )
 *
 * \brief	Callback for the multi threaded execution
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \param [in,out]	arg	If non-null, the argument.
 *
 * \return	.
 **************************************************************************************************/
template< class TJobData >
ITK_THREAD_RETURN_TYPE crlMultiThreadedWorkPileAbstract<TJobData>::WorkPileThreaderCallback( void *arg )
{
	//---------------------------------------
	// Get the parameters
	//---------------------------------------
	typedef itk::MultiThreader::ThreadInfoStruct  ThreadInfoType;
	ThreadInfoType * infoStruct = static_cast< ThreadInfoType * >( arg );
	Self *instance = (Self *)(infoStruct->UserData);
	const unsigned int threadId = infoStruct->ThreadID;

	//---------------------------------------
	// SYNCHRONIZATION OF ALL THE THREADS
	// (Wait that they all start)
	//---------------------------------------
	instance->m_StartBarrierSync->Wait();

	//---------------------------------------
	// Work on the workpile
	//---------------------------------------
	try {
		//---------------------------------------
		// Iterate over all the jobs
		//---------------------------------------
		TJobData newJob;
		while ( instance->GetNextJob(newJob) ) {
			instance->ExecuteJob(newJob, threadId);

			if ( instance->GetShowProgress() && instance->GetUseAutomaticProgressCounter() ) {
				instance->UpdateProgress(threadId);
			}
		}
	}
	catch (itk::ExceptionObject& e)
	{

		std::cout<< "THREAD ID"<<threadId<<" : ITK EXCEPTION ERROR CAUGHT"<<std::endl<< e.GetDescription() << std::endl << "Cannot continue." << std::endl ;
		throw e;
	}
	catch ( ... )
	{
		std::cout<<"THREAD ID"<<threadId<<" : UNKNOWN EXCEPTION ERROR." << std::endl << "Cannot continue."<< std::endl;
		throw;
	}

	//---------------------------------------
	// SYNCHRONIZATION OF ALL THE THREADS
	// (Wait that they all finish)
	//---------------------------------------
	instance->m_FinishBarrierSync->Wait();

	//---------------------------------------
	// Exit!
	//---------------------------------------
	return ITK_THREAD_RETURN_VALUE;
}


/**********************************************************************************************//**
 * \fn	template< class TInputImage, class TOutputImage > void WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::InitializeProgress(unsigned long totalCounter)
 *
 * \brief	Initializes the progress report.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \param	totalCounter	If not null, define the total counter.
 *
 * \return	.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::InitializeProgress(unsigned long totalCounter) 
{
	this->m_ProgressTimeProbe = itk::TimeProbe();
	this->m_CurrentProgress = 0;
	this->m_LastTimeProgressUpdated = 0;
	if ( totalCounter!=0 ) this->m_TotalProgress = totalCounter;
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPileAbstract<TJobData>::SetUseAutomaticProgressCounter(bool b)
 *
 * \brief	Sets whether the automatic progress counter is used (default: true). By default, the
 * 			progress counter is updated at each job. SetUseAutomaticProgressCounter(false) can be
 * 			used to update the counter more frequently (but manually) in order to show more
 * 			user-friendly remaining time.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 * \param	b	true to b.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::SetUseAutomaticProgressCounter(bool b)
{
	this->m_UseAutomaticProgressCounter = b;
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > bool crlMultiThreadedWorkPileAbstract<TJobData>::GetUseAutomaticProgressCounter() const
 *
 * \brief	Gets whether the automatic progress counter is used (default: true).By default, the
 * 			progress counter is updated at each job. SetUseAutomaticProgressCounter(false) can be
 * 			used to update the counter more frequently (but manually) in order to show more user-
 * 			friendly remaining time.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 *
 * \return	true if it succeeds, false if it fails.
 **************************************************************************************************/
template< class TJobData >
bool crlMultiThreadedWorkPileAbstract<TJobData>::GetUseAutomaticProgressCounter() const
{
	return this->m_UseAutomaticProgressCounter;
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPileAbstract<TJobData>::StartProgressTimeProbe(void)
 *
 * \brief	Starts the time probe used for showing remaining time during progress.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::StartProgressTimeProbe(void) 
{
	this->m_ProgressTimeProbe.Start();
	this->m_LastTimeProgressUpdated=0;
}

/**********************************************************************************************//**
 * \fn	template< class TInputImage, class TOutputImage > void WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::UpdateProgress( unsigned int threadId)
 *
 * \brief	Updates the progress report
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \param	threadId	Identifier for the thread.
 *
 * \return	.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPileAbstract<TJobData>::UpdateProgress( unsigned int threadId) 
{
	this->m_ProgressMutex.Lock();
	if ( this->m_CurrentProgress < this->m_TotalProgress )
	{
		this->m_CurrentProgress++;

		//-------------------------------------
		// Update cout not less than once per second
		//-------------------------------------
		m_ProgressTimeProbe.Stop();
		itk::TimeProbe::TimeStampType time = this->m_ProgressTimeProbe.GetTotal();
		m_ProgressTimeProbe.Start();

		if ( time-this->m_LastTimeProgressUpdated>1 )
		{
			this->m_LastTimeProgressUpdated = time;

#ifdef WIN32
			//-------------------------------------
			// If we are under win32, we can use 
			//functions to get/set the cursor pos
			//-------------------------------------
			CONSOLE_SCREEN_BUFFER_INFO csbi;
			COORD pos;
			if ( GetConsoleScreenBufferInfo( GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
				pos.Y = csbi.dwCursorPosition.Y-1;
			else
				pos.Y = 0;
			pos.X=0;
			SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pos);
#else
			//-------------------------------------
			// Linux/OSX: save current cursor position
			//-------------------------------------
			std::cout<<vt100_cursor_left(40); //"\e[s";
			std::cout<<vt100_cursor_up(1);
			//std::cout<< m_CurrentProgress << "/"  << m_TotalProgress <<std::endl;
#endif

			float percentage ;
			if ( this->m_TotalProgress==0 ) percentage=0.0f;
			else percentage = (float) ( this->m_CurrentProgress * 100.0 / ((double) this->m_TotalProgress));

			char szBuffer[128];
			sprintf(szBuffer,"%03d", (int)percentage);

			if ( this->m_CurrentProgress <= this->m_TotalProgress ) 
			{
				//cout<< d->m_LinePrefix.c_str() << szBuffer <<"%  "; 
				std::cout << "  " << szBuffer <<"%  "; 

				//-------------------------------------
				// If the task is finished, show computation time
				//-------------------------------------
				if ( this->m_CurrentProgress==this->m_TotalProgress ) 
				{
					this->m_ProgressTimeProbe.Stop();
					int t = (int)m_ProgressTimeProbe.GetTotal();
					std::string strTime = ConvertTimeToString(t);
					std::cout<< "("<<strTime <<")                  "; //<<  m_CurrentProgress << "/" << m_TotalProgress ;
				}
				//-------------------------------------	
				// Else show the estimated left time if needed and if possible (p>=1)
				//-------------------------------------
				else 
				{
					if ( percentage!=0.0f ) 
					{
						int t = (int)(100.0*time/percentage - time);

						if ( m_TimeFilter.size()>7 ) m_TimeFilter.pop_back();
						m_TimeFilter.insert(m_TimeFilter.begin(), t);

						double filteredTime=0;
						for ( unsigned int i=0; i<m_TimeFilter.size(); i++ )
							filteredTime+=m_TimeFilter[i];
						filteredTime=filteredTime/m_TimeFilter.size();

						if ( time<3 )
						{
							std::cout<< "(Estimating time left...)                ";
						}
						else
						{
							std::string strTime = ConvertTimeToString(filteredTime);
							std::cout<< "(Estimated left: "<<strTime.c_str()<<")                ";
						}
					}
				}

				std::cout << std::endl;
				std::flush(std::cout);
			}

		}
	}
	
	this->m_ProgressMutex.Unlock();

}

/**********************************************************************************************//**
 * \fn	template< class TInputImage, class TOutputImage > std::string WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::ConvertTimeToString(itk::RealTimeClock::TimeStampType t)
 *
 * \brief	Convert a time value (number of second) to a hh::mm::ss string.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \param	t	The.
 *
 * \return	.
 **************************************************************************************************/
template< class TJobData >
std::string crlMultiThreadedWorkPileAbstract<TJobData>::ConvertTimeToString(itk::RealTimeClock::TimeStampType t) 
{
	char szBuffer[512];

	int ss = ((int)t)%60; 
	int mm = (int)((t-ss)/60)%60;
	int hh = (int)(t / 3600);
	sprintf(szBuffer, "%02dh %02dm %02ds", hh, mm, ss );

	return std::string(szBuffer);
}



}

#endif
