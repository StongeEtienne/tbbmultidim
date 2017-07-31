#ifndef __itkTBBNDImageToImageFilter_h
#define __itkTBBNDImageToImageFilter_h

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    // We need that for TBB to correctly detect the Windows version
    //_WIN32_WINNT_WIN7
    //_WIN32_WINNT_WINXP
    //#define WINVER _WIN32_WINNT_WIN7
    //#define _WIN32_WINNT _WIN32_WINNT_WIN7
    #include "windows.h"
#endif

#include "itkImageToImageFilter.h"
#include "tbb/blocked_range.h"

/************************************************************************************************
* \class	TBBNDImageToImageFilter
*
* \brief	TBBNDImageToImageFilter.
*
* \author   Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
* \date     September 2016
*************************************************************************************************/

namespace itk{
	
	/**********************************************************************************************//**
	 * \class	TBBNDImageToImageFilter
	 *
	 * \brief   ImageToImageFilter using Intel Threading Building Blocks (TBB) parallelization
	 *          Multithreading with Thread and Job pool
	 *
	 * \author  Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	class  TBBNDFunctor;

	template< typename TInputImage, typename TOutputImage >
	class TBBNDImageToImageFilter : public ImageToImageFilter< TInputImage, TOutputImage>
	{
		friend class TBBNDFunctor<TInputImage,TOutputImage>;

	public:
		// Standard class typedefs.
		typedef TBBNDImageToImageFilter							Self;
		typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
		typedef SmartPointer< Self >							Pointer;
		typedef SmartPointer< const Self >						ConstPointer;

		// Run-time type information (and related methods).
		itkTypeMacro(TBBNDImageToImageFilter, ImageToImageFilter);

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

        unsigned int GetNbReduceDimensions() const;
		void SetNbReduceDimensions(int);

		virtual const ThreadIdType & GetNumberOfThreads() const;
        virtual void SetNumberOfThreads(ThreadIdType);

        // use TBBGenerateData instead of ThreadedGenerateData
        virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId) ITK_FINAL;
        virtual void TBBGenerateData(const OutputImageRegionType& outputRegionForThread); // todo pure virtual itkTypeMacro error

	protected:
		TBBNDImageToImageFilter();
		~TBBNDImageToImageFilter();

        void SetNumberOfJobs(unsigned int);
        unsigned int GetNumberOfJobs() const;

        void GenerateNumberOfJobs();

	private:
		TBBNDImageToImageFilter(const Self &); //purposely not implemented
		void operator=(const Self &); //purposely not implemented

		unsigned int m_TBBNumberOfJobs;
        unsigned int m_TBBNumberOfThreads;
        int m_TBBNbReduceDimensions;

        //struct ThreadStruct{ Pointer Filter;};
	};
}   //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTBBNDImageToImageFilter.txx"
#endif

#endif // __itkTBBNDImageToImageFilter_h
