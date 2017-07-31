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
#ifndef __crlMultiThreadedWorkPile_txx
#define __crlMultiThreadedWorkPile_txx

#include "crlMultiThreadedWorkPile.h"
#include <math.h>


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
crlMultiThreadedWorkPile<TJobData>::crlMultiThreadedWorkPile()
{

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
crlMultiThreadedWorkPile<TJobData>::~crlMultiThreadedWorkPile()
{

}


/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPile<TJobData>::AddJob(const TJobData& job);
 *
 * \brief	Adds a new job to the queue.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the work pile job.
 * \param	job	The job.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPile<TJobData>::QueueJob(const TJobData& job)
{
	this->LockJobQueueAccess();
	m_JobQueue.push_back(job);
	this->UnlockJobQueueAccess();
}

/**********************************************************************************************//**
 * \fn	template< class TJobData > void crlMultiThreadedWorkPile<TJobData>::ClearQueue()
 *
 * \brief	Clears the queue.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TJobData	Type of the job data.
 **************************************************************************************************/
template< class TJobData >
void crlMultiThreadedWorkPile<TJobData>::ClearQueue()
{
	this->LockJobQueueAccess();
	this->m_JobQueue.clear();
	this->UnlockJobQueueAccess();
}

/**********************************************************************************************//**
 * \fn	template< class TInputImage, class TOutputImage > bool WindowedMultiThreadedImageFilter< TInputImage, TOutputImage > ::GetNextRegionToProcess(OutputImageRegionType& nextRegion)
 *
 * \brief	Gets the next region to process.
 *
 * \author	Benoit Scherrer
 * \date	September 2012
 *
 * \param [in,out]	nextRegion	Reference to the next region.
 *
 * \return	True if there is another region to process. In that case modify nextRegion
 **************************************************************************************************/
template< class TJobData >
bool crlMultiThreadedWorkPile<TJobData>::GetNextJob( TJobData& job)
{
	bool success=true;

	this->LockJobQueueAccess();
	if ( m_JobQueue.size()==0 ) success=false;
	else {
		job = this->m_JobQueue.front();
		this->m_JobQueue.pop_front();
	}
	this->UnlockJobQueueAccess();

	return success;
}

// doc herited
template< class TJobData >
void crlMultiThreadedWorkPile<TJobData>::ExecuteQueue()
{
	int sizeQueue = 0;

	this->LockJobQueueAccess();
	sizeQueue = m_JobQueue.size();
	this->UnlockJobQueueAccess();

	if ( this->GetShowProgress() && this->GetUseAutomaticProgressCounter() ) 
		this->InitializeProgress(sizeQueue);
	
	Parent::ExecuteQueue();
}

}

#endif
