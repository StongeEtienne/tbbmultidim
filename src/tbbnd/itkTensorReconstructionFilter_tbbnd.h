#ifndef __TENSORRECONSTRUCTIONFILTER_TBBND_H
#define __TENSORRECONSTRUCTIONFILTER_TBBND_H

#include "itkTBBNDImageToImageFilter.h"
//#include "itkImageToImageFilter.h"

#include "itkDiffusionTensor3D.h"
#include "itkSpatialObject.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include <vector>
#include <Eigen/Dense>
#include "tbb/tick_count.h"

/************************************************************************************************
 * \class	TensorReconstruction
 *
 * \brief	TensorReconstruction.
 *
 * \author	Amir Jaberzadeh
 * \date	June 2015
*************************************************************************************************/

namespace itk
{
template< class TReferenceImagePixelType,
          class TGradientImagePixelType = TReferenceImagePixelType,
          class TTensorPixelType = double,
          class TMaskImageType = Image<unsigned char, 3> >                        
class ITK_EXPORT TensorReconstructionFilter_tbbnd:
        public TBBNDImageToImageFilter< Image< TReferenceImagePixelType, 3 >,
        Image< DiffusionTensor3D< TTensorPixelType >, 3 > > 
{
public:
    typedef TensorReconstructionFilter_tbbnd Self;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;
    typedef TBBNDImageToImageFilter< Image< TReferenceImagePixelType, 3 >,
    Image< DiffusionTensor3D< TTensorPixelType >, 3 >  >
    Superclass;
      
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    /** Runtime information support. */
    itkTypeMacro(TensorReconstructionFilter_tbbnd,
                 Superclass);
    typedef TReferenceImagePixelType ReferencePixelType;
    typedef TGradientImagePixelType GradientPixelType;
    typedef DiffusionTensor3D< TTensorPixelType > TensorPixelType;
    /** Reference image data, This image is acquired in the absence
* of a diffusion sensitizing field gradient */
    typedef typename Superclass::InputImageType ReferenceImageType;
    typedef Image< TensorPixelType, 3 > TensorImageType;
    typedef TensorImageType OutputImageType;
    typedef Image<TTensorPixelType, 3> OutputImageType2;
    typedef Image<TTensorPixelType, 3> OutputImageType3;
    typedef Image<TTensorPixelType, 3> OutputImageType4;

    typedef typename Superclass::OutputImageRegionType
    OutputImageRegionType;
    /** Typedef defining one (of the many) gradient images. */
    typedef Image< GradientPixelType, 3 > GradientImageType;
    /** An alternative typedef defining one (of the many) gradient images.
* It will be assumed that the vectorImage has the same dimension as the
* Reference image and a vector length parameter of \c n (number of
* gradient directions) */
    typedef VectorImage< GradientPixelType, 3 > GradientImagesType;

    /** The type for the optional SpatialObject for masking*/
    typedef SpatialObject<3> MaskSpatialObjectType;
    /** The type for the optional mask image */
    typedef TMaskImageType MaskImageType;
    /** Holds the tensor basis coefficients G_k */
    typedef Eigen::MatrixXd TensorBasisMatrixType;
    typedef Eigen::MatrixXd CoefficientMatrixType;
    /** Holds each magnetic field gradient used to acquire one DWImage */
    typedef vnl_vector_fixed< double, 3 > GradientDirectionType;
    /** Container to hold gradient directions of the 'n' DW measurements */
    typedef VectorContainer< TGradientImagePixelType,
    GradientDirectionType > GradientDirectionContainerType;

    typedef std::vector< double > VectorType;
    typedef enum {
        LLS,
        WLLS
    } solutionMethod;
    /** Set method to add a gradient direction and its corresponding image. */
    void AddGradientImage(const GradientDirectionType &, const GradientImageType *image);
    const GradientImageType *GetGradientImage(unsigned index) const;
    /** Another set method to add a gradient directions and its corresponding
* image. The image here is a VectorImage. The user is expected to pass the
* gradient directions in a container. The ith element of the container
* corresponds to the gradient direction of the ith component image the
* VectorImage. For the baseline image, a vector of all zeros
* should be set. */
    void SetGradientImage(GradientDirectionContainerType *,
                          const GradientImagesType *image);
    /** Set method to set the reference image. */
    void SetReferenceImage(ReferenceImageType *referenceImage)
    {
        if ( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
        {
            itkExceptionMacro(<< "Cannot call both methods:"
                              << "AddGradientImage and SetGradientImage. Please call only one of them.");
        }
        this->ProcessObject::SetNthInput(0, referenceImage);
        m_GradientImageTypeEnumeration = GradientIsInManyImages;
    }
    /** Get reference image */
    virtual ReferenceImageType * GetReferenceImage()
    { return ( static_cast< ReferenceImageType * >( this->ProcessObject::GetInput(0) ) ); }
    /** Return the gradient direction. idx is 0 based */
    virtual GradientDirectionType GetGradientDirection(TGradientImagePixelType idx) const
    {
        if ( idx >= m_NumberOfGradientDirections )
        {
            itkExceptionMacro(<< "Gradient direction " << idx << "does not exist");
        }
        return m_GradientDirectionContainer->ElementAt(idx + 1);
    }
    /** set an image mask */
    void SetMaskImage(MaskImageType *maskImage);
    /** set a spatial object mask */
    void SetMaskSpatialObject(MaskSpatialObjectType *maskSpatialObject);
    /** Threshold on the reference image data. The output tensor will be a null
* tensor for pixels in the reference image that have a value less than this
* threshold. */
    itkSetMacro(Threshold, ReferencePixelType);
    itkGetConstMacro(Threshold, ReferencePixelType);
    /**
* The BValue \f$ (s/mm^2) \f$ value used in normalizing the tensors to
* physically meaningful units. See equation (24) of the first reference for
* a description of how this is applied to the tensor estimation.
* Equation (1) of the same reference describes the physical significance.
*/
     void SetBValue( const std::vector<double>& v) { m_BValue = v;}

#ifdef GetBValue
#undef GetBValue
#endif
      itkGetMacro(BValue, std::vector<double>);
//    itkGetConstReferenceMacro(BValue, VectorType);
    /** set the solution method */
    itkSetMacro( Method, solutionMethod );
#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro( ReferenceEqualityComparableCheck,
                     ( Concept::EqualityComparable< ReferencePixelType > ) );
    itkConceptMacro( TensorEqualityComparableCheck,
                     ( Concept::EqualityComparable< TensorPixelType > ) );
    itkConceptMacro( GradientConvertibleToDoubleCheck,
                     ( Concept::Convertible< GradientPixelType, double > ) );
    itkConceptMacro( DoubleConvertibleToTensorCheck,
                     ( Concept::Convertible< double, TensorPixelType > ) );
    itkConceptMacro( GradientReferenceAdditiveOperatorsCheck,
                     ( Concept::AdditiveOperators< GradientPixelType, GradientPixelType,
                       ReferencePixelType > ) );
    itkConceptMacro( GradientReferenceAdditiveAndAssignOperatorsCheck,
                     ( Concept::AdditiveAndAssignOperators< GradientPixelType,
                       ReferencePixelType > ) );
    itkConceptMacro( ReferenceOStreamWritableCheck,
                     ( Concept::OStreamWritable< ReferencePixelType > ) );
    itkConceptMacro( TensorOStreamWritableCheck,
                     ( Concept::OStreamWritable< TensorPixelType > ) );
    /** End concept checking */
#endif
    OutputImageType* GetOutput1(){
        return dynamic_cast< OutputImageType * >(this->ProcessObject::GetOutput(0) );
    }
    OutputImageType2* GetOutput2(){
        return dynamic_cast< OutputImageType2 * >(this->ProcessObject::GetOutput(1) );
    }
    OutputImageType3* GetOutput3(){
        return dynamic_cast< OutputImageType3 * >(this->ProcessObject::GetOutput(2) );
    }
    OutputImageType4* GetOutput4(){
        return dynamic_cast< OutputImageType4 * >(this->ProcessObject::GetOutput(3) );
    }

    void setResidualImage(){
        m_residImage = true;
    }
    void setADCImage(){
        m_ADCImage = true;
    }
    double m_beforetime;
protected:
    TensorReconstructionFilter_tbbnd();
    ~TensorReconstructionFilter_tbbnd() {}
    void PrintSelf(std::ostream & os, Indent indent) const;
    void ComputeTensorBasis();
    void BeforeThreadedGenerateData();
    void TBBGenerateData(const
                              OutputImageRegionType & outputRegionForThread);
    void AfterThreadedGenerateData();
    /** enum to indicate if the gradient image is specified as a single multi-
* component image or as several separate images */
    typedef enum {
        GradientIsInASingleImage = 1,
        GradientIsInManyImages,
        Else
    } GradientImageTypeEnumeration;
private:
    /* Tensor basis coeffs */
    TensorBasisMatrixType m_TensorBasis;
    CoefficientMatrixType m_BMatrix;
    /** container to hold gradient directions */
    typename GradientDirectionContainerType::Pointer m_GradientDirectionContainer;
    /** Number of gradient measurements */
    unsigned int m_NumberOfGradientDirections;
    /** Number of baseline images */
    unsigned int m_NumberOfBaselineImages;
    /** Number of all images */
    unsigned int m_NumImages;	
    /** Threshold on the reference image data */
    ReferencePixelType m_Threshold;
    /** LeBihan's b-value for normalizing tensors */
    std::vector<double> m_BValue;
    std::vector<double> newBValue;
    /** Gradient image was specified in a single image or in multiple images */
    GradientImageTypeEnumeration m_GradientImageTypeEnumeration;
    /** Mask Image Present */
    bool m_MaskImagePresent;
    solutionMethod m_Method;
    bool m_residImage;
    bool m_ADCImage;
};
} // end of namespace itk


#endif // __TENSORRECONSTRUCTIONFILTER_TBBND_H

