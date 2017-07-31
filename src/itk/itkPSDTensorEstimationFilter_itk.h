#ifndef __ITKPSDTENSORESTIMATIONFILTER_ITK_H
#define __ITKPSDTENSORESTIMATIONFILTER_ITK_H

//#include "itkTBBImageToImageFilter.h"
#include "itkImageToImageFilter.h"

#include "itkSpatialObject.h"
#include "itkDiffusionTensor3D.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"
#include <math.h>
#include <vector>
#include <list>

#include "itkSingleValuedCostFunction.h"

/************************************************************************************************
 * \class	PSDTensorEstimationFilter_itk
 *
 * \brief	PSDTensorEstimationFilter_itk
 *
 * \author	Amir Jaberzadeh, Etienne St-Onge
 * \date	June 2017
*************************************************************************************************/

namespace itk
{
template< class TInputImagePixelType,
          class TOutputImageType,
          class TTensorPixelType = double,
          class TMaskImageType = Image<unsigned char, 3> >
class ITK_EXPORT PSDTensorEstimationFilter_itk:
        public ImageToImageFilter< Image< TInputImagePixelType, 3 >,
        Image< DiffusionTensor3D< TTensorPixelType >, 3 > >
{
public:
    typedef PSDTensorEstimationFilter_itk Self;
    typedef SmartPointer< Self > Pointer;
    typedef SmartPointer< const Self > ConstPointer;
    typedef ImageToImageFilter< Image< TInputImagePixelType, 3 >,
    Image< DiffusionTensor3D< TTensorPixelType >, 3 > >
    Superclass;
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    /** Runtime information support. */
    itkTypeMacro(PSDTensorEstimationFilter_itk,
                 Superclass);
    typedef TInputImagePixelType ReferencePixelType;
    typedef TInputImagePixelType GradientPixelType;
    typedef DiffusionTensor3D< TTensorPixelType > TensorPixelType;
    /** Reference image data, This image is acquired in the absence
* of a diffusion sensitizing field gradient */
    typedef typename Superclass::InputImageType ReferenceImageType;
    typedef Image< TensorPixelType, 3 > TensorImageType;
    typedef TensorImageType OutputImageType1;
    typedef Image<float, 3> OutputImageType2;
    typedef Image<double, 3> OutputImageType3;
    typedef Image< TInputImagePixelType, 3 > OutputImageType4;
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
    typedef vnl_matrix_fixed< double, 6, 6 > TensorBasisMatrixType;
    typedef vnl_matrix< double > CoefficientMatrixType;
    /** Holds each magnetic field gradient used to acquire one DWImage */
    typedef vnl_vector_fixed< double, 3 > GradientDirectionType;
    /** Container to hold gradient directions of the 'n' DW measurements */
    typedef VectorContainer< TInputImagePixelType,
    GradientDirectionType > GradientDirectionContainerType;
    /** Set method to add a gradient direction and its corresponding image. */
    typedef enum {
          CLLS,
          CWLLS,
          CNLS,
          CWNLS,
          LTS
      } solutionMethod;
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
    virtual GradientDirectionType GetGradientDirection(unsigned int idx) const
    {
        if ( idx >= m_NumberOfGradientDirections )
        {
            itkExceptionMacro(<< "Gradient direction " << idx << "does not exist");
        }
        return m_GradientDirectionContainer->ElementAt(idx + 1);
    }
    /** set an image mask */
    void SetMaskImage(MaskImageType *maskImage);
    /** set an initial tensor image */
    void SetInitialTensorImage(const TensorImageType *InitialTensorImage){
        m_InitialTensor = const_cast<TensorImageType *>(InitialTensorImage);
        m_TensorImagePresent = true;
    }
    void SetInitialB0Image(const Image< TInputImagePixelType, 3 > *InitialB0Image){
    	m_B0Image = const_cast<Image< TInputImagePixelType, 3 > *>(InitialB0Image);
    	m_B0ImagePresent = true;
    }
    /** set a spatial object mask */
    void SetMaskSpatialObject(MaskSpatialObjectType *maskSpatialObject);
    /** Threshold on the reference image data. The output tensor will be a null
* tensor for pixels in the reference image that have a value less than this
* threshold. */
    itkSetMacro(Threshold, float);
    itkGetConstMacro(Threshold, float);
    /**
* The BValue \f$ (s/mm^2) \f$ value used in normalizing the tensors to
* physically meaningful units. See equation (24) of the first reference for
* a description of how this is applied to the tensor estimation.
* Equation (1) of the same reference describes the physical significance.
*/
    void SetBValue( const std::vector<double>& v) { m_BValue = v; }
#ifdef GetBValue
#undef GetBValue
#endif
     itkGetMacro(BValue, std::vector<double>);
    /** set the solution method */
    itkSetMacro( Method, solutionMethod );
    void SetMeanResiduals(){m_MeanResiduals = true;}
    void SetAllResiduals(){m_AllResiduals = true;}
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
    OutputImageType1* GetOutput1(){
        return dynamic_cast< OutputImageType1 * >(this->ProcessObject::GetOutput(0) );
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
    void SetMaxIterations(unsigned short Maxiter){
    m_Maxiter = Maxiter;
    }
   double m_beforetime;
protected:
    PSDTensorEstimationFilter_itk();
    ~PSDTensorEstimationFilter_itk() {}
//    DataObject::Pointer MakeOutput(unsigned int idx);
    void PrintSelf(std::ostream & os, Indent indent) const;
    void ComputeTensorBasis();
    void BeforeThreadedGenerateData();
    void AfterThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, 
                              			itk::ThreadIdType);
    /** enum to indicate if the gradient image is specified as a single multi-
* component image or as several separate images */
    typedef enum {
        GradientIsInASingleImage = 1,
        GradientIsInManyImages,
        Else
    } GradientImageTypeEnumeration;
    bool m_MeanResiduals;
    bool m_AllResiduals;
private:
    /* Tensor basis coeffs */
    TensorBasisMatrixType m_TensorBasis;
    TensorImageType* m_InitialTensor;
    Image< TInputImagePixelType, 3 > * m_B0Image;
    bool m_TensorImagePresent;
    CoefficientMatrixType m_BMatrix;
    /** container to hold gradient directions */
    typename GradientDirectionContainerType::Pointer m_GradientDirectionContainer;
    /** Number of gradient measurements */
    unsigned int m_NumberOfGradientDirections;
    /** Number of baseline images */
    unsigned int m_NumberOfBaselineImages;
    /** Threshold on the reference image data */
    float m_Threshold;
    /** LeBihan's b-value for normalizing tensors */
    std::vector<double> m_BValue;
    /** Gradient image was specified in a single image or in multiple images */
    GradientImageTypeEnumeration m_GradientImageTypeEnumeration;
    /** Mask Image Present */
    bool m_MaskImagePresent;
    solutionMethod m_Method;
    unsigned short m_Maxiter;
    unsigned int m_NumImages;
    bool m_B0ImagePresent;

//    typename OutputImageType3::Pointer residImage;
};

template<class TInputPixelType, class TOutputType>
class ITK_EXPORT MyCostFunction1
    : public itk::SingleValuedCostFunction
{
public:
  /** Standard class typedefs */
  typedef MyCostFunction1					Self;
  typedef SingleValuedCostFunction			Superclass;
  typedef SmartPointer< Self >				Pointer;
  typedef SmartPointer< const Self >		ConstPointer;
  typedef vnl_matrix< double > CoefficientMatrixType;
  typedef itk::PSDTensorEstimationFilter_itk<
            TInputPixelType, TOutputType >  PSDTensorEstimationFilterType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MyCostFunction1, SingleValuedCostFunction) ;

  /** Method for creation through the object factory. */
  itkNewMacro(Self) ;
  itkSetMacro( Method, typename PSDTensorEstimationFilterType::solutionMethod );
  virtual unsigned int GetNumberOfParameters() const { return 7; }
  void SetDesignMatrix(const CoefficientMatrixType* BMatrix){
        m_BMatrix = const_cast<CoefficientMatrixType*> (BMatrix);
  }

  void SetPixelLogValue(const vnl_vector< double >* B){
      m_PixelLog = const_cast<vnl_vector< double >* > (B);
  }
  void SetPixelIntensityValue(const std::vector< double >* B){
      m_PixelValue = const_cast<std::vector< double >* > (B);
  }

  /** This method returns the value of the cost function corresponding
    * to the specified parameters. */
  virtual MeasureType GetValue( const ParametersType & parameters ) const
  {
      MeasureType sum = 0;
      MeasureType inner_sum = 0;

      if(m_Method == PSDTensorEstimationFilterType::CLLS){
          for(int i = 0; i < m_PixelLog->size(); i++){
              inner_sum = m_BMatrix->get(0,i) * (parameters[0]*parameters[0]) +
                      m_BMatrix->get(1,i) * parameters[0]*parameters[3] +
                      m_BMatrix->get(2,i) * parameters[0]*parameters[5] +
                      m_BMatrix->get(3,i) * (parameters[1]*parameters[1] + parameters[3]*parameters[3]) +
                      m_BMatrix->get(4,i) * (parameters[1]*parameters[4] + parameters[3]*parameters[5]) +
                      m_BMatrix->get(5,i)* (parameters[2]*parameters[2] + parameters[4]*parameters[4] + parameters[5]*parameters[5]);
              sum = pow(m_PixelLog->get(i) -(vcl_log(parameters[6] * 1e5) + inner_sum), 2) + sum;     // cost func for CLLS method
          }
      }
      else if(m_Method == PSDTensorEstimationFilterType::CWLLS){
          for(int i = 0; i < m_PixelLog->size(); i++){
              inner_sum = m_BMatrix->get(0,i) * (parameters[0]*parameters[0]) +
                      m_BMatrix->get(1,i) * parameters[0]*parameters[3] +
                      m_BMatrix->get(2,i) * parameters[0]*parameters[5] +
                      m_BMatrix->get(3,i) * (parameters[1]*parameters[1] + parameters[3]*parameters[3]) +
                      m_BMatrix->get(4,i) * (parameters[1]*parameters[4] + parameters[3]*parameters[5]) +
                      m_BMatrix->get(5,i)* (parameters[2]*parameters[2] + parameters[4]*parameters[4] + parameters[5]*parameters[5]);
              sum = pow(m_PixelValue->at(i),2) * pow(m_PixelLog->get(i) - (vcl_log(parameters[6] * 1e5) + inner_sum), 2) + sum;   // cost func for CWLLS method
          }
      }
      else if(m_Method == PSDTensorEstimationFilterType::CNLS){
          for(int i = 0; i < m_PixelLog->size(); i++){
              inner_sum = m_BMatrix->get(0,i) * (parameters[0]*parameters[0]) +
                      m_BMatrix->get(1,i) * parameters[0]*parameters[3] +
                      m_BMatrix->get(2,i) * parameters[0]*parameters[5] +
                      m_BMatrix->get(3,i) * (parameters[1]*parameters[1] + parameters[3]*parameters[3]) +
                      m_BMatrix->get(4,i) * (parameters[1]*parameters[4] + parameters[3]*parameters[5]) +
                      m_BMatrix->get(5,i) * (parameters[2]*parameters[2] + parameters[4]*parameters[4] + parameters[5]*parameters[5]);
              sum = pow((m_PixelValue->at(i) - (parameters[6] * 1e5) * vcl_exp(inner_sum)), 2) + sum;     // cost func for CNLS method
          }
      }
      else if(m_Method == PSDTensorEstimationFilterType::CWNLS){
          for(int i = 0; i < m_PixelLog->size(); i++){
              inner_sum = m_BMatrix->get(0,i) * (parameters[0]*parameters[0]) +
                      m_BMatrix->get(1,i) * parameters[0]*parameters[3] +
                      m_BMatrix->get(2,i) * parameters[0]*parameters[5] +
                      m_BMatrix->get(3,i) * (parameters[1]*parameters[1] + parameters[3]*parameters[3]) +
                      m_BMatrix->get(4,i) * (parameters[1]*parameters[4] + parameters[3]*parameters[5]) +
                      m_BMatrix->get(5,i)* (parameters[2]*parameters[2] + parameters[4]*parameters[4] + parameters[5]*parameters[5]);
              sum = pow(m_PixelValue->at(i),2) * pow((m_PixelValue->at(i) - parameters[6] * 1e5 * vcl_exp(inner_sum)), 2) + sum;     // cost func for CWNLS method
          }
      }
      else if(m_Method == PSDTensorEstimationFilterType::LTS){
          std::list<MeasureType> mylist;
          std::list<MeasureType>::iterator it;
          mylist.clear();
          for(int i = 0; i < m_PixelLog->size(); i++){
              inner_sum = m_BMatrix->get(0,i) * (parameters[0]*parameters[0]) +
                      m_BMatrix->get(1,i) * parameters[0]*parameters[3] +
                      m_BMatrix->get(2,i) * parameters[0]*parameters[5] +
                      m_BMatrix->get(3,i) * (parameters[1]*parameters[1] + parameters[3]*parameters[3]) +
                      m_BMatrix->get(4,i) * (parameters[1]*parameters[4] + parameters[3]*parameters[5]) +
                      m_BMatrix->get(5,i)* (parameters[2]*parameters[2] + parameters[4]*parameters[4] + parameters[5]*parameters[5]);
               mylist.push_back(pow((m_PixelValue->at(i) - vcl_log(parameters[6] * 1e5) * vcl_exp(inner_sum)), 2));
          }
          mylist.sort();
          it = mylist.begin();
          for(int i = 0; i < m_PixelLog->size()/1.5; i++, it++)  //To do: change trimming parameter 1.5
              sum = *it +sum;
//              std::cout<< *it <<" ";
      }
      return 0.5*sum;
  }

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.
    NOT USED IN NEWUOA */
  virtual void GetDerivative( const ParametersType &,
                                    DerivativeType & ) const
  {
  }


protected:
  MyCostFunction1() {

  }
  virtual ~MyCostFunction1() {}
  virtual void PrintSelf(std::ostream& os, Indent indent) const
  {}

private:
  CoefficientMatrixType* m_BMatrix;
  vnl_vector< double >* m_PixelLog;
  std::vector< double >* m_PixelValue;
  typename PSDTensorEstimationFilterType::solutionMethod m_Method;
}; // end of class
}  //end namespace itk

#endif // __ITKPSDTENSORESTIMATIONFILTER_H
