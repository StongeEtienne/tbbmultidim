#ifndef ITKRICIANNOISECORRECTIONFILTER_ITKND_H
#define ITKRICIANNOISECORRECTIONFILTER_ITKND_H

//#include "itkTBBImageToImageFilter.h"
#include "itkCRLImageToImageFilter.h"

#include <vector>

/************************************************************************************************
 * \class	RicianNoiseCorrectionFilter_itknd
 *
 * \brief	RicianNoiseCorrectionFilter_itknd
 *
 * \author	Amir Jaberzadeh, Etienne St-Onge
 * \date	June 2017
*************************************************************************************************/

namespace itk{
using namespace std;
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
class RicianNoiseCorrectionFilter_itknd : public CRLImageToImageFilter< TInputImage, TOutputImage>
{
public:
   typedef RicianNoiseCorrectionFilter_itknd Self;
   typedef CRLImageToImageFilter< TInputImage, TOutputImage > Superclass;
   typedef SmartPointer< Self > Pointer;
   typedef SmartPointer< const Self > ConstPointer;

   typedef TInputImage InputImageType;
   typedef typename InputImageType::Pointer InputImagePointer;

   typedef TOutputImage OutputImageType;
   typedef typename OutputImageType::Pointer OutputImagePointer;
   typedef typename OutputImageType::RegionType OutputImageRegionType;
   typedef typename OutputImageType::PixelType OutputImagePixelType;

   typedef TMaskImage MaskImageType;
   typedef typename MaskImageType::Pointer MaskImagePointer;
   typedef typename TMaskImage::PixelType MaskPixelType;

    itkNewMacro(Self);
    itkTypeMacro(RicianNoiseCorrectionFilter_itknd, Superclass);

   virtual void BeforeThreadedGenerateData();
   void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, 
                              			itk::ThreadIdType);
   void SetMaskImage(const TMaskImage *mask);
   void SetVariance(float variance){
   	m_Var = variance;
   	m_MaskPresent = false;
   }
   const TMaskImage * GetMaskImage() const;
   float CalculateMean(vector<float> imageVector);
   float CalculateVariance(vector<float> imageVector);
   double m_beforetime;

protected:
    RicianNoiseCorrectionFilter_itknd();
    ~RicianNoiseCorrectionFilter_itknd();
private:
    RicianNoiseCorrectionFilter_itknd(const Self &); //purposely not implemented
    void operator=(const Self &); //purposely not implemented
    std::vector<float> m_Variance;
    float m_Var;
    bool m_MaskPresent;
};
}  // end namespace itk

#endif // ITKRICIANNOISECORRECTIONFILTER_ITKND_H
