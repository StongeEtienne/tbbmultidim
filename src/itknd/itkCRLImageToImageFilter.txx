/*
 * Copyright (c) 2016 Benoit Scherrer, Etienne Saint-Onge, Boston Children's Hospital.
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
*/

#ifndef __itkCRLImageToImageFilter_hxx
#define __itkCRLImageToImageFilter_hxx

#include "itkCRLImageToImageFilter.h"
#include "itkImageSource.h"
#include "itkImageRegionSplitterBase.h"
#include "itkOutputDataObjectIterator.h"
#include "vnl/vnl_math.h"

/************************************************************************************************
 * \class	CRLImageToImageFilter
 *
 * \brief	CRLImageToImageFilter.
 *
 * \author  Benoit Scherrer and Etienne St-Onge
 * \date	September 2016
*************************************************************************************************/

namespace itk {

	// Constructor
	template< typename TInputImage, typename TOutputImage >
	CRLImageToImageFilter< TInputImage, TOutputImage >::CRLImageToImageFilter()
	{
		// By default, Automatic NbReduceDimensions
        this->SetNbReduceDimensions(-1);
		this->m_StartBarrierSync = itk::Barrier::New();
		this->m_FinishBarrierSync = itk::Barrier::New();
	}

	// Destructor
	template< typename TInputImage, typename TOutputImage >
	CRLImageToImageFilter< TInputImage, TOutputImage >::~CRLImageToImageFilter()
	{
	}



	/************************************************************************************************
	 * \fn	template< typename TInputImage, typename TOutputImage > void CRLImageToImageFilter< TInputImage, TOutputImage > ::GenerateData()
	 *
	 * \brief	New default implementation for GenerateData() to use TBB
	 *
	 * \author  Benoit Scherrer and Etienne St-Onge
	 *
	 * \tparam	typename TInputImage 	Type of the typename t input image.
	 * \tparam	typename TOutputImage	Type of the typename t output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	void CRLImageToImageFilter< TInputImage, TOutputImage >::GenerateData()
	{
        // Get the size of the requested region
        typename TOutputImage::ConstPointer output = static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
        typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();

        // Call a method that can be overriden by a subclass to allocate
        // memory for the filter's outputs
        this->AllocateOutputs();

        // Generate the number of Jobs
        // based on the OutputImageDimension, NumberOfThreads and NbReduceDimensions
        this->GenerateNumberOfJobs();
		
		// Reinitialize current job index
		this->ResetJobQueue();

		// Call a method that can be overridden by a subclass to perform
        // some calculations prior to splitting the main computations into
        // separate threads
		this->BeforeThreadedGenerateData();

		// Initializes the barrier
		this->m_StartBarrierSync->Initialize(this->GetNumberOfThreads());
		this->m_FinishBarrierSync->Initialize(this->GetNumberOfThreads());

		// Set up the multithreaded processing
		this->GetMultiThreader()->SetNumberOfThreads( this->GetNumberOfThreads() );
		this->GetMultiThreader()->SetSingleMethod(this->MyThreaderCallback, (void *)this );

		// multithread the execution
		this->GetMultiThreader()->SingleMethodExecute();

		// Call a method that can be overridden by a subclass to perform
		// some calculations after all the threads have completed
		this->AfterThreadedGenerateData();
	}

	/**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > ITK_THREAD_RETURN_TYPE CRLImageToImageFilter< TInputImage, TOutputImage >::MyThreaderCallback( void *arg )
	 *
	 * \brief	Internal function. Callback method for the multithreader.
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2016
	 *
	 * \exception	Thrown an exception if an error occurs.
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 * \param [in,out]	Poiner to the itk::MultiThreader::ThreadInfoStruct
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	ITK_THREAD_RETURN_TYPE CRLImageToImageFilter< TInputImage, TOutputImage >::MyThreaderCallback( void *arg )
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
		int jobId;
		try {
			while ( (jobId=instance->GetNextJob())>=0 ) 
			{
				instance->ExecuteJob(jobId);
			}
		}
		catch (itk::ExceptionObject& e)
		{

			std::cout<< "THREAD ID"<<threadId<<" / JOB ID << " << jobId << ": ITK EXCEPTION ERROR CAUGHT"<<std::endl<< e.GetDescription() << std::endl << "Cannot continue." << std::endl ;
			throw e;
		}
		catch ( ... )
		{
			std::cout<<"THREAD ID"<<threadId<<" / JOB ID << " << jobId << " : UNKNOWN EXCEPTION ERROR." << std::endl << "Cannot continue."<< std::endl;
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
	 * \fn	template< typename TInputImage, typename TOutputImage > int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNextJob()
	 *
	 * \brief	Gets the next job.
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2016
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 *
	 * \return	The next job&lt;typename t input image,typename t output image &gt;
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNextJob()
	{
		int jobId = -1;

		this->m_JobQueueMutex.Lock();
		if ( m_CurrentJobQueueIndex == m_TBBNumberOfJobs ) jobId=-1;
		else 
		{
			jobId = m_CurrentJobQueueIndex;
			m_CurrentJobQueueIndex++;
		}
		this->m_JobQueueMutex.Unlock();

		return jobId;
	}

	/**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > int CRLImageToImageFilter< TInputImage, TOutputImage >::ExecuteJob( int jobId )
	 *
	 * \brief	Executes the job identified by jobId.
	 *
	 * \author	Benoit Scherrer
	 * \date	October 2016
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 * \param	jobId	Identifier for the job.
	 *
	 * \return	.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	void CRLImageToImageFilter< TInputImage, TOutputImage >::ExecuteJob( int jobId )
	{
		// Get the size of the requested region
		typename TOutputImage::ConstPointer output = static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
		typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();
		typename TOutputImage::SizeType size = outputSize;

		//typename TOutputImage::SizeType size = m_OutputSize;

		typename TOutputImage::IndexType index;
		index.Fill(0);

		if (this->GetNbReduceDimensions() > 0)
		{
			unsigned int i = OutputImageDimension - (unsigned int)this->GetNbReduceDimensions();

			index[i] = jobId;
			size[i] = 1;
			while (i < OutputImageDimension - 1)
			{
				index[i+1] = index[i] / outputSize[i];
				index[i] = index[i] % outputSize[i];
				size[i+1] = 1;
				i++;
			}
		}

		// Construct an itk::ImageRegion
		OutputImageRegionType myRegion(index, size);

		// Run the ThreadedGenerateData method! 
		this->ThreadedGenerateData(myRegion, jobId);
	}

	/************************************************************************************************
	 * \fn	template< typename TInputImage, typename TOutputImage > unsigned int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
	 *
	 * \brief	Gets the number of jobs.
	 * 			
	 * 			\warning This function only returns a valid value when called from AllocateOutputs,
	 * 			BeforeThreadedGenerateData(), TBBGenerateData() and AfterThreadedGenerateData(),
	 * 			ie when the input image is known
	 *
	 * \author  Benoit Scherrer and Etienne St-Onge
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 *
	 * \return	The number of jobs
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	unsigned int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
	{
		return m_TBBNumberOfJobs;
	}


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > unsigned int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNbReduceDimensions() const
     *
     * \brief   Gets the number of dimension to separate for the Jobs multithreading
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	unsigned int CRLImageToImageFilter< TInputImage, TOutputImage >::GetNbReduceDimensions() const
	{
		return m_TBBNbReduceDimensions;
	}


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void CRLImageToImageFilter< TInputImage, TOutputImage >::SetNbReduceDimensions(int nbReduceDim)
     *
     * \brief   Set the number of dimension to separate and multithread each section.
     *          (nbReduceDim <= 0  : negative number for automatic splitting)
     *
     *          \example : for a 3D image (volume) with the shape 30x10x5
     *          nbReduceDim == 0  : Will generate a single (1) Job with the whole image (size 30x10x5)
     *          nbReduceDim == 1  : Will generate 5 Jobs with the slices (size 30x10)
     *          nbReduceDim == 2  : Will generate 50 Jobs with the lines (size 30)
     *          nbReduceDim == 3  : Will generate 1500 Jobs with each voxel (size 1)
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	void CRLImageToImageFilter< TInputImage, TOutputImage >::SetNbReduceDimensions(int nbReduceDim)
	{
        if (nbReduceDim > (int)OutputImageDimension)
        {
            this->m_TBBNbReduceDimensions = (int)OutputImageDimension;
        }
        else
        {
            this->m_TBBNbReduceDimensions = nbReduceDim;
        }
	}


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void CRLImageToImageFilter< TInputImage, TOutputImage >::GenerateNumberOfJobs()
     *
     * \brief   Generate the number Jobs based on the NbReduceDimensions
     *              or based on the NumberOfThreads and the Image Dimension (if NbReduceDimensions was not set).
     *
     *          \warning  This function must be called after the NumberOfThreads is set.
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void CRLImageToImageFilter< TInputImage, TOutputImage >::GenerateNumberOfJobs()
    {
        // Get the size of the requested region
        typename TOutputImage::ConstPointer output = static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
        typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();

        // Generate the number of job
        if (m_TBBNbReduceDimensions < 0)
        {
            // assert (GetNumberOfThreads()>0)
            // This function must be called after the NumberOfThreads is Set

            // Heuristic NbReduceDimensions
            m_TBBNbReduceDimensions = 0;
            m_TBBNumberOfJobs = 1;
            int current_dim = OutputImageDimension-1;

            // Minimum Number of Jobs, based on the Number of thread
            unsigned int MinNbJobs = 20*this->GetNumberOfThreads();
            while( current_dim >= 0 && m_TBBNumberOfJobs < MinNbJobs )
            {
                ++m_TBBNbReduceDimensions;
                m_TBBNumberOfJobs *= outputSize[current_dim];
                --current_dim;
            }
        }
        else
        {
            // Fixed (preset NbReduceDimensions)
            m_TBBNumberOfJobs = 1;
            for (unsigned int i = OutputImageDimension - GetNbReduceDimensions(); i < OutputImageDimension; ++i)
            {
                m_TBBNumberOfJobs *= outputSize[i];
            }
        }
    }

   template< typename TInputImage, typename TOutputImage >
    void CRLImageToImageFilter< TInputImage, TOutputImage >::ResetJobQueue()
    {
		this->m_CurrentJobQueueIndex = 0;
	}

    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void CRLImageToImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
     *
     * \brief   Use *TBBGenerateData()* instead of ThreadedGenerateData with CRLImageToImageFilter
     *
     *          \warning CRLImageToImageFilter doesn't support threadId
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    //template< typename TInputImage, typename TOutputImage >
    //void CRLImageToImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
    //{
    //    std::cout << " ERROR : Subclass should override this method!!! (redefine *ThreadedGenerateData*)" << std::endl;
    //}



}  //namespace itk

#endif
