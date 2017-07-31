#ifndef ITKRICIANNOISECORRECTIONFILTER_ITKND_TXX
#define ITKRICIANNOISECORRECTIONFILTER_ITKND_TXX

#include "itkRicianNoiseCorrectionFilter_itknd.h"
#include "itkImageRegionIterator.h"
#include "itkCRLImageToImageFilter.txx"

#include <math.h>
#include "tbb/tick_count.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/************************************************************************************************
 * \class	RicianNoiseCorrectionFilter_itknd
 *
 * \brief	RicianNoiseCorrectionFilter_itknd
 *
 * \author	Amir Jaberzadeh, Etienne St-Onge
 * \date	June 2017
*************************************************************************************************/

namespace itk {


template< typename TInputImage, typename TMaskImage,typename TOutputImage >
void
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType &outputRegionForThread,
itk::ThreadIdType ThreadId)
{
    //     cout<<"this was the threaded version "<<endl;

    InputImagePointer inputPtr = const_cast< TInputImage * >( this->GetInput() );
    OutputImagePointer outputPtr = this->GetOutput();
    itk::ImageRegionIterator<InputImageType> inputIterator(inputPtr, outputRegionForThread);
    itk::ImageRegionIterator<OutputImageType> outputIterator(outputPtr, outputRegionForThread);

    inputIterator.GoToBegin();
    outputIterator.GoToBegin();

    std::vector<float> background_noise_stddev;
    //     std::cout<<" m_variance "<<m_Variance.size()<<std::endl;
    for(int i = 0; i < m_Variance.size();i++)
        background_noise_stddev.push_back(sqrt(m_Variance[i]));
    while(!inputIterator.IsAtEnd())
    {
        typename InputImageType::PixelType pixval;
        typename InputImageType::PixelType outval = inputIterator.Get();
        for(int i = 0; i < m_Variance.size(); i++){
            pixval = inputIterator.Get();
            //         for (unsigned int i = 0; i < 1000; i++) {  // Used to increase granularity
            if( pixval[i] > background_noise_stddev[i]) {
                outval[i] = sqrt( (pixval[i]*pixval[i]) - m_Variance[i] );
            } else {
                outval[i] = 0.0;
            }
            outputIterator.Set(outval);
            //         }
        }
        ++inputIterator;
        ++outputIterator;
    }
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::SetMaskImage(const TMaskImage *mask)
{
    this->ProcessObject::SetNthInput( 1, const_cast< TMaskImage * >( mask ) );
}
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
const TMaskImage *
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::GetMaskImage() const
{
    return itkDynamicCastInDebugMode< MaskImageType * >( const_cast< DataObject * >( this->ProcessObject::GetInput(1) ) );
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::RicianNoiseCorrectionFilter_itknd()
{
	m_MaskPresent = true;
	m_Var = 0;

}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::~RicianNoiseCorrectionFilter_itknd()
{

}


template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::BeforeThreadedGenerateData()
{
    tbb::tick_count starttime, endtime; // for timing
    starttime = tbb::tick_count::now();
    typename InputImageType::ConstPointer input = this->GetInput(0);    
    itk::ImageRegionConstIterator<InputImageType> InputIterator( input, input->GetLargestPossibleRegion() );
    InputIterator.GoToBegin();
    typename InputImageType::PixelType pixel = InputIterator.Get();
    if(m_MaskPresent){
    MaskImagePointer maskPtr =
            const_cast< TMaskImage * >( this->GetMaskImage() );
    itk::ImageRegionConstIterator<MaskImageType> MaskIterator( maskPtr, maskPtr->GetLargestPossibleRegion() );
    MaskIterator.GoToBegin();

    float result;vector<float> ImageVector;
    for(int i = 0; i < pixel.GetSize(); i++){
        InputIterator.GoToBegin();
        MaskIterator.GoToBegin();
        while(!MaskIterator.IsAtEnd())
        {
            // Set the current pixel to white
            pixel = InputIterator.Get();
            result = MaskIterator.Get() * pixel[i];
            if(result!=0)
                ImageVector.push_back(result);
            ++MaskIterator;
            ++InputIterator;
        }
        m_Variance.push_back(this->CalculateVariance(ImageVector));
        ImageVector.clear();
    }
    }
    else{
    	for(int i = 0; i < pixel.GetSize(); i++)
    		m_Variance.push_back(m_Var);
    }
    endtime = tbb::tick_count::now();
    m_beforetime = (endtime - starttime).seconds();
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
float
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::CalculateMean(vector<float> imageVector)
{
    int max = imageVector.size();
    //    cout<<"max is "<<max<<" ";
    float sum = 0;
    for(int i = 0; i < max; i++){
        sum += imageVector[i];
        //        cout<<sum<<" ";
    }
    return  (sum / (float) max);
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
float
RicianNoiseCorrectionFilter_itknd< TInputImage, TMaskImage, TOutputImage >
::CalculateVariance(vector<float> imageVector)
{
    float mean = this->CalculateMean(imageVector);

        //  cout<< "mean is "<<mean<<endl;
    float temp = 0.0;
    int max = imageVector.size();
    for(int i = 0; i < max; i++)
    {
        temp += (imageVector[i] - mean) * (imageVector[i] - mean) ;
    }
    //      cout<<temp<<" ";
    float variance = (temp / max)/(2.0 - M_PI/2.0); // Gudbjartsson 1995, The Rician Distribution of Noisy MRI Data eq 3
    // m_Variance = 0;
   //     cout << " variance is " << variance<<std::endl;
    return variance;
}
}  // end namespace itk

#endif // ITKRICIANNOISECORRECTIONFILTER_ITK_TXX
