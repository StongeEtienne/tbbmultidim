
#ifndef __crlDilateLabelImageFilter2_txx
#define __crlDilateLabelImageFilter2_txx

#include "crlDilateLabelImageFilter2.h"
#include <limits.h>

#include "itkConstantBoundaryCondition.h"
#include "itkNumericTraits.h"
#include "itkMorphologyImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"

template<class TInputImage, class TOutputImage, class TKernel>
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::crlDilateLabelImageFilter2()
{
  m_DefaultBoundaryCondition.SetConstant( NumericTraits<PixelType>::NonpositiveMin() );
  m_BoundaryCondition = &m_DefaultBoundaryCondition;
}


template <class TInputImage, class TOutputImage, class TKernel>
void 
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer  inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  
  if ( !inputPtr )
    {
    return;
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Kernel.GetRadius() );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template<class TInputImage, class TOutputImage, class TKernel>
void
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::BeforeThreadedGenerateData()
{
	Superclass::BeforeThreadedGenerateData();

	this->m_DistanceImage = DistanceImageType::New();
	this->m_DistanceImage->SetOrigin(this->GetInput()->GetOrigin() );
	this->m_DistanceImage->SetSpacing(this->GetInput()->GetSpacing() );
	this->m_DistanceImage->SetDirection(this->GetInput()->GetDirection() );
	this->m_DistanceImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
	this->m_DistanceImage->Allocate();

	this->m_DistanceImage->FillBuffer(9999999);
}



template<class TInputImage, class TOutputImage, class TKernel>
void
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId) 
{
  // Neighborhood iterators
  NeighborhoodIteratorType b_iter;

  // Find the boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> fC;
  faceList = fC(this->GetInput(), outputRegionForThread, m_Kernel.GetRadius());

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;

  ImageRegionIterator<TOutputImage> o_iter;
  DistanceImageIteratorType d_iter;

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Process the boundary faces, these are N-d regions which border the
  // edge of the buffer

  const KernelIteratorType kernelBegin = m_Kernel.Begin();
  const KernelIteratorType kernelEnd = m_Kernel.End();
  
  for (fit = faceList.begin(); fit != faceList.end(); ++fit)
    { 
    b_iter = NeighborhoodIteratorType(m_Kernel.GetRadius(),
                                      this->GetInput(), *fit);
    
    o_iter = ImageRegionIterator<OutputImageType>(this->GetOutput(), *fit);
	d_iter = ImageRegionIterator<DistanceImageType>(this->m_DistanceImage, *fit);

	b_iter.OverrideBoundaryCondition(m_BoundaryCondition);
    b_iter.GoToBegin();

    while ( ! o_iter.IsAtEnd() )
      {
      o_iter.Set( this->Evaluate(b_iter, kernelBegin, kernelEnd, d_iter) );
      ++b_iter;
      ++o_iter;
      progress.CompletedPixel();
      }
    }
  
}

template<class TInputImage, class TOutputImage, class TKernel>
void
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Kernel: " << m_Kernel << std::endl;
  os << indent << "Boundary condition: " << typeid( *m_BoundaryCondition ).name() << std::endl;
}



template<class TInputImage, class TOutputImage, class TKernel>
typename crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>::PixelType
crlDilateLabelImageFilter2<TInputImage, TOutputImage, TKernel>
::Evaluate(const NeighborhoodIteratorType &nit,
           const KernelIteratorType kernelBegin,
           const KernelIteratorType kernelEnd,
		   DistanceImageIteratorType& distImageIt)
{
  unsigned int i;
  float dst;
  PixelType max = NumericTraits<PixelType>::NonpositiveMin();

  typename InputImageType::IndexType centerIndex, neighbIndex;

  KernelIteratorType kernel_it;
	
  PixelType nCenPxVal=nit.GetCenterPixel();
	if ( nCenPxVal!=0 ) return nCenPxVal;

  centerIndex = nit.GetIndex();

  for( i=0, kernel_it=kernelBegin; kernel_it<kernelEnd; ++kernel_it, ++i )
    {
    // if structuring element is positive, use the pixel under that element
    // in the image
    if( *kernel_it > NumericTraits<KernelPixelType>::Zero )
      {
      // note we use GetPixel() on the SmartNeighborhoodIterator to
      // respect boundary conditions
      temp = nit.GetPixel(i);
	 

	  neighbIndex = nit.GetIndex (i);
	  dst=0;
	  for ( unsigned k=0; k<ImageDimension; k++ )
		  dst+=(neighbIndex[k]-centerIndex[k])*(neighbIndex[k]-centerIndex[k]);

	  if ( dst < distImageIt.Value() )
	  {
		  
	  }
		if( temp  max )
        {
        max = temp;
        }
      }
    }
  
  return max;
} 


#endif

