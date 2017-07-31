#ifndef __crlDilateLabelImageFilter2_h
#define __crlDilateLabelImageFilter2_h

#include "itkImageToImageFilter.h"

#include "itkNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"
#include "itkConstSliceIterator.h"
#include "itkImageBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "itkSmartPointer.h"


using namespace itk;

/**********************************************************************************************//**
 * \class	crlDilateLabelImageFilter2
 *
 * \brief	Crl label dilate image filter. 
 *
 * \author	Benoit Scherrer
 * \date	February 2011
*************************************************************************************************/
template<class TInputImage, class TOutputImage, class TKernel>
class crlDilateLabelImageFilter2 : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard Self typedef */
  typedef crlDilateLabelImageFilter2                    Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;
  
  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(crlDilateLabelImageFilter2, ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef TInputImage                                   InputImageType;
  typedef TOutputImage                                  OutputImageType;
  typedef typename TInputImage::RegionType              RegionType;
  typedef typename TInputImage::SizeType                SizeType;
  typedef typename TInputImage::IndexType               IndexType;
  typedef typename TInputImage::PixelType               PixelType;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  
  typedef typename itk::Image<float, TInputImage::ImageDimension> DistanceImageType;
  typedef itk::ImageRegionIterator<DistanceImageType> DistanceImageIteratorType;

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Typedef for boundary conditions. */
  typedef ImageBoundaryCondition<InputImageType> *       ImageBoundaryConditionPointerType;
  typedef ImageBoundaryCondition<InputImageType> const * ImageBoundaryConditionConstPointerType;
  typedef ConstantBoundaryCondition<InputImageType>      DefaultBoundaryConditionType;
  

/** Neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<TInputImage> 
  NeighborhoodIteratorType;

  /** Kernel typedef. */
  typedef TKernel KernelType;
  
  /** Kernel (structuring element) iterator. */
  typedef typename KernelType::ConstIterator KernelIteratorType;
  
  /** n-dimensional Kernel radius. */
  typedef typename KernelType::SizeType RadiusType;

  /** Set kernel (structuring element). */
  itkSetMacro(Kernel, KernelType);

  /** Get the kernel (structuring element). */
  itkGetConstReferenceMacro(Kernel, KernelType);
  
  /** Type of the pixels in the Kernel. */
  typedef typename TKernel::PixelType            KernelPixelType;

  /** MorphologyImageFilters need to make sure they request enough of an
   * input image to account for the structuring element size.  The input
   * requested region is expanded by the radius of the structuring element.
   * If the request extends past the LargestPossibleRegion for the input,
   * the request is cropped by the LargestPossibleRegion. */
  void GenerateInputRequestedRegion();

  /** Allows a user to override the internal boundary condition. Care should be
   * be taken to ensure that the overriding boundary condition is a persistent
   * object during the time it is referenced.  The overriding condition
   * can be of a different type than the default type as long as it is
   * a subclass of ImageBoundaryCondition. */
  void OverrideBoundaryCondition(const ImageBoundaryConditionPointerType i)
    { m_BoundaryCondition = i; }

  /** Rest the boundary condition to the default */
  void ResetBoundaryCondition()
    { m_BoundaryCondition = &m_DefaultBoundaryCondition; }
  
  /** Get the current boundary condition. */
  itkGetConstMacro(BoundaryCondition, ImageBoundaryConditionPointerType);
  
protected:
  crlDilateLabelImageFilter2();
  ~crlDilateLabelImageFilter2() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData();

  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData (const OutputImageRegionType& 
                              outputRegionForThread,
                              int threadId);

  /** Evaluate image neighborhood with kernel to find the new value 
   * for the center pixel value. */
  virtual PixelType Evaluate(const NeighborhoodIteratorType &nit,
                             const KernelIteratorType kernelBegin,
                             const KernelIteratorType kernelEnd,
							 DistanceImageIteratorType& distImageIt);

private:
  crlDilateLabelImageFilter2(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** kernel or structuring element to use. */
  KernelType m_Kernel;

  /** Pointer to a persistent boundary condition object used
   * for the image iterator. */
  ImageBoundaryConditionPointerType m_BoundaryCondition;

  /** Default boundary condition */
  DefaultBoundaryConditionType m_DefaultBoundaryCondition;
  
  typename DistanceImageType::Pointer m_DistanceImage;
}; // end of class


#ifndef ITK_MANUAL_INSTANTIATION
#include "crlDilateLabelImageFilter2.txx"
#endif

#endif

