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

#ifndef __itkCRLImageToImageFilter_h
#define __itkCRLImageToImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkBarrier.h"

/************************************************************************************************
* \class	CRLImageToImageFilter
*
* \brief	CRLImageToImageFilter.
*
* \author   Benoit Scherrer and Etienne St-Onge
* \date     September 2016
*************************************************************************************************/

namespace itk{
	
	/**********************************************************************************************//**
	 * \class	CRLImageToImageFilter
	 *
	 * \brief   ImageToImageFilter using a job pool
	 *
	 * \author Benoit Scherrer and Etienne St-Onge
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	class CRLImageToImageFilter : public ImageToImageFilter< TInputImage, TOutputImage>
	{
	public:
		// Standard class typedefs.
		typedef CRLImageToImageFilter							Self;
		typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
		typedef SmartPointer< Self >							Pointer;
		typedef SmartPointer< const Self >						ConstPointer;

		// Run-time type information (and related methods).
		itkTypeMacro(CRLImageToImageFilter, ImageToImageFilter);

		// Superclass typedefs.
		typedef typename Superclass::OutputImageRegionType	OutputImageRegionType;
		typedef typename Superclass::OutputImagePixelType	OutputImagePixelType;

		// Some convenient typedefs. The same as itk::ImageToImageFilter
		typedef TInputImage									InputImageType;
		typedef typename InputImageType::Pointer			InputImagePointer;
		typedef typename InputImageType::ConstPointer		InputImageConstPointer;
		typedef typename InputImageType::RegionType			InputImageRegionType;
		typedef typename InputImageType::PixelType			InputImagePixelType;
		typedef typename TInputImage::SizeType				InputImageSizeType;

		typedef TOutputImage								OutputImageType;
		typedef typename TOutputImage::SizeType				OutputImageSizeType;


		// ImageDimension constants
		itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
		itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

	public:
		itkNewMacro(Self);

		void GenerateData();

        unsigned int	GetNbReduceDimensions() const;
		void			SetNbReduceDimensions(int);

	protected:
		CRLImageToImageFilter();
		~CRLImageToImageFilter();

        unsigned int	GetNumberOfJobs() const;
        void			GenerateNumberOfJobs();

		int				GetNextJob();
		void			ExecuteJob( int jobId );
		void			ResetJobQueue() ;

		static ITK_THREAD_RETURN_TYPE MyThreaderCallback( void *arg );

	protected:
		typename itk::Barrier::Pointer			m_StartBarrierSync;
		typename itk::Barrier::Pointer			m_FinishBarrierSync;

	private:
		CRLImageToImageFilter(const Self &); //purposely not implemented
		void operator=(const Self &); //purposely not implemented

		int										m_CurrentJobQueueIndex;
		itk::SimpleFastMutexLock				m_JobQueueMutex;

		unsigned int	m_TBBNumberOfJobs;
        unsigned int	m_TBBNumberOfThreads;
        int				m_TBBNbReduceDimensions;

        //struct ThreadStruct{ Pointer Filter;};
	};
}   //end namespace itk


#endif // __itkCRLImageToImageFilter_h
