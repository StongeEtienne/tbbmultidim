#ifndef __itkTBBImageToImageFilter_hxx
#define __itkTBBImageToImageFilter_hxx


#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
  // We need that for TBB to correctly detect the Windows version
  #include "windows.h"
#endif
#include "tbb/tbb_stddef.h"
#include "itkTBBImageToImageFilter.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/compat/thread"
#include <sstream>
/************************************************************************************************
 * \class	TBBImageToImageFilter
 *
 * \brief	TBBImageToImageFilter.
 *
 * \author	Amir Jaberzadeh and Benoit Scherrer
 * \date	May 2015
*************************************************************************************************/

namespace itk {

	/**********************************************************************************************//**
	 * \class	TbbFunctor
	 *
	 * \brief	Tbb functor to execute jobs in parallel.
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	class TBBFunctor 
	{
	public:
		typedef TBBFunctor		Self;
		typedef TOutputImage	OutputImageType;
		typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
		typedef typename TOutputImage::SizeType OutputImageSizeType;
		typedef typename OutputImageType::RegionType OutputImageRegionType;

		itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
		itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

		typedef TBBImageToImageFilter<TInputImage,TOutputImage> TbbImageFilterType;

		TBBFunctor(TbbImageFilterType *tbbFilter, const OutputImageSizeType& outputSize):
		m_TbbFilter(tbbFilter), m_OutputSize(outputSize) {}

		void operator() ( const tbb::blocked_range<int>& r ) const
		{
			// Setup the size of the jobs to be done
			typename TOutputImage::SizeType size = m_OutputSize;
			size[OutputImageDimension - 1] = r.end() - r.begin();

			// Setup the starting index
			typename TOutputImage::IndexType index;
			index.Fill(0);
			index[OutputImageDimension - 1] = r.begin();

			// Construct an itk::ImageRegion
			OutputImageRegionType myRegion(index, size);

			// Run the ThreadedGenerateData method!
			m_TbbFilter->TBBGenerateData(myRegion);

		}

	private:
		TbbImageFilterType *m_TbbFilter;
		OutputImageSizeType m_OutputSize;
	};



	// Constructor
	template< typename TInputImage, typename TOutputImage >
	TBBImageToImageFilter< TInputImage, TOutputImage >::TBBImageToImageFilter()
	{
		// By default, do not define the number of threads.
		// Let TBB doing that.
		this->SetNumberOfThreads(0);
		m_NumberOfThreads = false;
	}

	// Destructor
	template< typename TInputImage, typename TOutputImage >
	TBBImageToImageFilter< TInputImage, TOutputImage >::~TBBImageToImageFilter()
	{ }

	/**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > void TBBImageToImageFilter< TInputImage, TOutputImage > ::GenerateData()
	 *
	 * \brief	New default implementation for GenerateData() to use TBB
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	typename TInputImage 	Type of the typename t input image.
	 * \tparam	typename TOutputImage	Type of the typename t output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	void TBBImageToImageFilter< TInputImage, TOutputImage >::GenerateData()
	{
		// Get the size of the requested region
		typename TOutputImage::ConstPointer output = static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
		typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();
		this->m_NumberOfJobs = outputSize[OutputImageDimension - 1];

		// Call a method that can be overriden by a subclass to allocate
		// memory for the filter's outputs
		this->AllocateOutputs();

		// Call a method that can be overridden by a subclass to perform
		// some calculations prior to splitting the main computations into
		// separate threads
		this->BeforeThreadedGenerateData();

		// Set up the number of threads. Only for testing purposes. Should not
		// be used in practice.
		tbb::task_scheduler_init init(-2);
		if (m_NumberOfThreads)
			init.initialize(this->GetNumberOfThreads());
		else
			init.initialize();
		
		// Do the task decomposition using parallel_for
		tbb::parallel_for( tbb::blocked_range<int>(0, this->m_NumberOfJobs, 1),
			TBBFunctor<TInputImage, TOutputImage>(this, outputSize));

		// Call a method that can be overridden by a subclass to perform
		// some calculations after all the threads have completed
		this->AfterThreadedGenerateData();
	}
	
	 /**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > void TBBImageToImageFilter< TInputImage, TOutputImage > ::ThreadedGenerateData()
	 *
	 * \brief	finalize usage of ThreadedGenerateData() to use TBBGenerateData() instead
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	typename TInputImage 	Type of the typename t input image.
	 * \tparam	typename TOutputImage	Type of the typename t output image.
	 **************************************************************************************************/
	 template< typename TInputImage, typename TOutputImage > 
	 void TBBImageToImageFilter< TInputImage, TOutputImage > 
	 ::ThreadedGenerateData(const OutputImageRegionType&, itk::ThreadIdType){
	// The following code is equivalent to:
	// itkExceptionMacro("subclass should override this method!!!");
	// The ExceptionMacro is not used because gcc warns that a
	// 'noreturn' function does return
 	std::ostringstream message;
	message << "itk::ERROR: " << this->GetNameOfClass()
	<< "(" << this << "): " << "Subclass should override TBBGenerateData() method instead!!!";
	ExceptionObject e_(__FILE__, __LINE__, message.str().c_str(),ITK_LOCATION);
	throw e_;
	 }
	 /**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > void TBBImageToImageFilter< TInputImage, TOutputImage > ::TBBGenerateData()
	 *
	 * \brief	Overload implementation of TBBGenerateData()
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	typename TInputImage 	Type of the typename t input image.
	 * \tparam	typename TOutputImage	Type of the typename t output image.
	 **************************************************************************************************/
	 template< typename TInputImage, typename TOutputImage > 
	 void TBBImageToImageFilter< TInputImage, TOutputImage > 
	 ::TBBGenerateData(const OutputImageRegionType&){
	// The following code is equivalent to:
	// itkExceptionMacro("subclass should override this method!!!");
	// The ExceptionMacro is not used because gcc warns that a
	// 'noreturn' function does return
 	std::ostringstream message;
	message << "itk::ERROR: " << this->GetNameOfClass()
	<< "(" << this << "): " << "Subclass should override this method!!!";
	ExceptionObject e_(__FILE__, __LINE__, message.str().c_str(),ITK_LOCATION);
	throw e_;
	}
	 	 
	/**********************************************************************************************//**
	 * \fn	template< typename TInputImage, typename TOutputImage > unsigned int TBBImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
	 *
	 * \brief	Gets the number of jobs.
	 * 			
	 * 			\warning This function only returns a valid value when called from AllocateOutputs,
	 * 			BeforeThreadedGenerateData, ThreadedGenerateData and AfterThreadedGenerateData,
	 * 			ie when the input image is known
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 *
	 * \return	The number of jobs
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	unsigned int TBBImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
	{
		return m_NumberOfJobs;
	}

}  //namespace itk

#endif
