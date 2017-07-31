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

#ifndef __crlMultiThreadedWorkPileAbstract_h
#define __crlMultiThreadedWorkPileAbstract_h

#include <itkMultiThreader.h>
#include <itkBarrier.h>
#include <itkSimpleFastMutexLock.h>
#include <itkTimeProbe.h>
#include <vector>

namespace crl {


#ifndef myGetSetMacros
#define myGetSetMacros(name,type) \
	virtual void Set##name (const type _arg) \
	{ \
	this->m_##name = _arg; \
    } \
	virtual type Get##name () \
	{ \
	return this->m_##name; \
    }
#endif

/**********************************************************************************************//**
 * \class	crlMultiThreadedWorkPileAbstract
 *
 * \brief	Define a threaded workpile
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 **************************************************************************************************/
template< class TJobData>
class crlMultiThreadedWorkPileAbstract
{
public:
	typedef crlMultiThreadedWorkPileAbstract		Self;
	typedef TJobData								JobData;

	/** Set/Get if show progress  */
	myGetSetMacros( ShowProgress, bool );
	myGetSetMacros( NbThreads, int );

	itk::MultiThreader * GetMultiThreader() const;

	crlMultiThreadedWorkPileAbstract();
	~crlMultiThreadedWorkPileAbstract();

	// Virtual function to inmplement
	virtual bool GetNextJob(TJobData& job) = 0;

	virtual void ExecuteQueue();

	void	LockJobQueueAccess();
	void	UnlockJobQueueAccess();

protected:
	// Trick to provide a pointer to the current object to the
	// thread call-back function, that can only be static 
	// (and cannot use 'this')
	/*struct ThreadedWorkPileDataStruct
	{
		Self::Pointer	instance;
	};*/

	void		InitializeProgress(unsigned long totalCounter=0) ;
	void		SetUseAutomaticProgressCounter(bool b);
	bool		GetUseAutomaticProgressCounter() const;
	void		StartProgressTimeProbe(void) ;
	void		UpdateProgress(unsigned int threadId) ;
	std::string	ConvertTimeToString(itk::RealTimeClock::TimeStampType t) ;

	static ITK_THREAD_RETURN_TYPE WorkPileThreaderCallback( void *arg );

	virtual void ExecuteJob( const TJobData &job, int threadId) = 0;

protected:
	typename itk::MultiThreader::Pointer	m_MultiThreader;
	typename itk::Barrier::Pointer			m_StartBarrierSync;
	typename itk::Barrier::Pointer			m_FinishBarrierSync;

private:
	bool									m_ShowProgress;

	int										m_NbThreads;

	itk::SimpleFastMutexLock				m_JobQueueMutex;

	bool									m_UseAutomaticProgressCounter;
	unsigned long							m_CurrentProgress;
	unsigned long							m_TotalProgress;
	itk::TimeProbe							m_ProgressTimeProbe;
	itk::SimpleFastMutexLock				m_ProgressMutex;
	int										m_LastTimeProgressUpdated;

	std::vector<int>						m_TimeFilter;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlMultiThreadedWorkPileAbstract.txx"
#endif

#endif
