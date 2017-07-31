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

#ifndef __crlMultiThreadedWorkPile_h
#define __crlMultiThreadedWorkPile_h

#include "crlMultiThreadedWorkPileAbstract.h"
#include <deque>

namespace crl {

/**********************************************************************************************//**
 * \class	crlThreadedWorkPile
 *
 * \brief	Define a threaded workpile
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 **************************************************************************************************/
template< class TJobData>
class crlMultiThreadedWorkPile : public crlMultiThreadedWorkPileAbstract<TJobData>
{
public:
	typedef crlMultiThreadedWorkPile					Self;
	typedef crlMultiThreadedWorkPileAbstract<TJobData>	Parent;
	typedef TJobData									JobData;


	crlMultiThreadedWorkPile();
	~crlMultiThreadedWorkPile();

	void	ClearQueue();
	void	QueueJob(const TJobData& job);
	virtual bool GetNextJob(TJobData& job);

	virtual void ExecuteQueue();

protected:
	// to implement: virtual void ExecuteJob( const TJobData &job, int threadId) = NULL;
	
protected:
	std::deque<TJobData>				m_JobQueue;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "crlMultiThreadedWorkPile.txx"
#endif

#endif
