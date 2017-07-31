#ifndef __itkTBBImageToImageFilter_h
#define __itkTBBImageToImageFilter_h

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
  // We need that for TBB to correctly detect the Windows version
  #include "windows.h"
#endif

#include "itkImageToImageFilter.h"
#include "tbb/blocked_range.h"

/************************************************************************************************
* \class	TBBImageToImageFilter
*
* \brief	TBBImageToImageFilter.
*
* \author	Amir Jaberzadeh and Benoit Scherrer
* \date	May 2015
*************************************************************************************************/

namespace itk{
	
	/**********************************************************************************************//**
	 * \class	TBBImageToImageFilter
	 *
	 * \brief	New itk::TBBImageToImageFilter
	 *
	 * \author	Amir Jaberzadeh and Benoit Scherrer
	 *
	 * \tparam	TInputImage 	Type of the input image.
	 * \tparam	TOutputImage	Type of the output image.
	 **************************************************************************************************/
	template< typename TInputImage, typename TOutputImage >
	class  TBBFunctor;

	template< typename TInputImage, typename TOutputImage >
	class TBBImageToImageFilter : public ImageToImageFilter< TInputImage, TOutputImage>
	{
		friend class TBBFunctor<TInputImage,TOutputImage>;

	public:
		/** Standard class typedefs. */
		typedef TBBImageToImageFilter							Self;
		typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
		typedef SmartPointer< Self >							Pointer;
		typedef SmartPointer< const Self >						ConstPointer;

		/** Run-time type information (and related methods). */
		itkTypeMacro(TBBImageToImageFilter, ImageToImageFilter);

		/** Superclass typedefs. */
		typedef typename Superclass::OutputImageRegionType	OutputImageRegionType;
		typedef typename Superclass::OutputImagePixelType	OutputImagePixelType;

		/** Some convenient typedefs. The same as itk::ImageToImageFilter */
		typedef TInputImage									InputImageType;
		typedef typename InputImageType::Pointer			InputImagePointer;
		typedef typename InputImageType::ConstPointer		InputImageConstPointer;
		typedef typename InputImageType::RegionType			InputImageRegionType;
		typedef typename InputImageType::PixelType			InputImagePixelType;
		typedef typename TInputImage::SizeType				InputImageSizeType;

		/** ImageDimension constants */
		itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
		itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

	public:
		itkNewMacro(Self);
		void GenerateData();
		void SetManualThreads(){
			m_NumberOfThreads = true;
		}
		
	protected:
		TBBImageToImageFilter();
		~TBBImageToImageFilter();
		virtual void ThreadedGenerateData(const OutputImageRegionType&, itk::ThreadIdType);
		virtual void TBBGenerateData(const OutputImageRegionType&);
		unsigned int GetNumberOfJobs() const;
		
	private:
		TBBImageToImageFilter(const Self &); //purposely not implemented
		void operator=(const Self &); //purposely not implemented

		unsigned int m_NumberOfJobs;
		bool m_NumberOfThreads;
	};
}   //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTBBImageToImageFilter.txx"
#endif

#endif // __itkTBBImageToImageFilter_h
