#ifndef __TENSORRECONSTRUCTIONFILTER_TBBND_TXX
#define __TENSORRECONSTRUCTIONFILTER_TBBND_TXX

#include "itkTensorReconstructionFilter_tbbnd.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkArray.h"
#include "itkImageMaskSpatialObject.h"
#include "vnl/vnl_vector.h"
#include <math.h>
#include "itkImageFileWriter.h"
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
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::TensorReconstructionFilter_tbbnd()
{
    // At least 1 inputs is necessary for a vector image.
    // For images added one at a time we need at least six
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(3);
    this->SetNthOutput( 0, ( OutputImageType::New() ).GetPointer());
    this->SetNthOutput( 1, ( OutputImageType2::New() ).GetPointer());
    this->SetNthOutput( 2, ( OutputImageType3::New() ).GetPointer());
    this->SetNthOutput( 3, ( OutputImageType4::New() ).GetPointer());
    m_NumberOfGradientDirections = 0;
    m_NumberOfBaselineImages = 1;
    m_NumImages = 1;
    m_Threshold = NumericTraits< ReferencePixelType >::min();
    m_GradientImageTypeEnumeration = Else;
    m_GradientDirectionContainer = NULL;
    m_TensorBasis.setIdentity();
//    m_BValue = NULL;
    m_MaskImagePresent = false;
    m_Method = LLS;
    m_residImage = false;
    m_ADCImage = false;
}


template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::BeforeThreadedGenerateData()
{
    tbb::tick_count starttime, endtime; // for timing
    starttime = tbb::tick_count::now();

    // If we have more than 3 inputs, then each input, except the first is a
    // gradient image. The number of gradient images must match the number of
    // gradient directions.
    const unsigned int numberOfInputs = this->GetNumberOfIndexedInputs();
    // There need to be at least 6 gradient directions to be able to compute the
    // tensor basis
    if ( m_NumberOfGradientDirections < 6 )
    {
        itkExceptionMacro(<< "At least 6 gradient directions are required");
    }

    // If there is only 1 gradient image, it must be an itk::VectorImage.
    // Otherwise
    // we must have a container of (numberOfInputs-1) itk::Image. Check to make
    // sure
    if ( (numberOfInputs == 1 ||
          (numberOfInputs == 2 && this->m_MaskImagePresent))
         && m_GradientImageTypeEnumeration != GradientIsInASingleImage )
    {
        std::string gradientImageClassName(
                    this->ProcessObject::GetInput(0)->GetNameOfClass() );
        if ( strcmp(gradientImageClassName.c_str(), "VectorImage") != 0 )
        {
            itkExceptionMacro(
                        << "There is only one Gradient image. I expect that to be a VectorImage. "
                        << "But its of type: " << gradientImageClassName);
        }
    }
    this->ComputeTensorBasis();

    endtime = tbb::tick_count::now();
    m_beforetime = (endtime - starttime).seconds();

    // If there's a mask, make sure it matches the dimensions of the
    // gradient image(s).
    if(!this->m_MaskImagePresent)
    {
        return;
    }
    typename ImageMaskSpatialObject<3>::Pointer maskSpatialObject =
            dynamic_cast<ImageMaskSpatialObject<3> *>(this->ProcessObject::GetInput(1));
    if(maskSpatialObject.IsNull())
    {
        return; // not a mask image
    }
    typename MaskImageType::ConstPointer maskImage = maskSpatialObject->GetImage();
    typename MaskImageType::SizeType maskSize =
            maskImage->GetLargestPossibleRegion().GetSize();
    typename MaskImageType::SizeType refSize;
    typename MaskImageType::PointType maskOrigin =
            maskImage->GetOrigin();
    typename MaskImageType::PointType refOrigin;
    typename MaskImageType::SpacingType maskSpacing =
            maskImage->GetSpacing();
    typename MaskImageType::SpacingType refSpacing;
    typename MaskImageType::DirectionType maskDirection =
            maskImage->GetDirection();
    typename MaskImageType::DirectionType refDirection;
    if( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        ReferenceImageType *refImage =
                static_cast<ReferenceImageType *>(this->ProcessObject::GetInput(0));
        refSize =
                refImage->GetLargestPossibleRegion().GetSize();
        refOrigin = refImage->GetOrigin();
        refSpacing = refImage->GetSpacing();
        refDirection = refImage->GetDirection();
    }
    else if( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
        GradientImagesType *gradientImage =
                static_cast< GradientImagesType * >(this->ProcessObject::GetInput(0) );
        refSize =
                gradientImage->GetLargestPossibleRegion().GetSize();
        refOrigin = gradientImage->GetOrigin();
        refSpacing = gradientImage->GetSpacing();
        refDirection = gradientImage->GetDirection();
    }
    // size mismatch is a deal breaker. Iterators are useless.
    if(refSize != maskSize)
    {
        itkExceptionMacro( << "Mask size doesn't match Reference Image Size"
                           << " Mask Size " << maskSize
                           << " Ref Size " << refSize );
    }
/*    
    // Origin, Spacing, Direction, should match but it isn't fatal if
    // they don't.
    
    if(refOrigin != maskOrigin)
    {
        itkWarningMacro( << "Mask origin doesn't match Reference origin "
                         << "Mask Origin " << maskOrigin
                         << " Ref Origin " << refOrigin );
    }
    if(refSpacing != maskSpacing)
    {
        itkWarningMacro( << "Mask spacing doesn't match Reference spacing "
                         << "Mask Spacing " << maskSpacing
                         << " Ref Spacing " << refSpacing );
    }
    if(refDirection != maskDirection)
    {
        itkWarningMacro( << "Mask direction doesn't match Reference direction "
                         << "Mask Direction " << maskDirection
                         << " Ref Direction " << refDirection );
    }
*/
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::TBBGenerateData(const OutputImageRegionType & outputRegionForThread)
{
    typename OutputImageType::Pointer outputImage =
            this->GetOutput1();
    typename OutputImageType2::Pointer outputImage2 =
            this->GetOutput2();
    typename OutputImageType3::Pointer outputImage3 =
            this->GetOutput3();
    typename OutputImageType4::Pointer outputImage4 =
            this->GetOutput4();
    ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
    ImageRegionIterator< OutputImageType2 > resiImageIt(outputImage2, outputRegionForThread);
    ImageRegionIterator< OutputImageType3 > ADCImageIt(outputImage3, outputRegionForThread);
    ImageRegionIterator< OutputImageType4 > b0ImageIt(outputImage4, outputRegionForThread);

    oit.GoToBegin();
    resiImageIt.GoToBegin();
    ADCImageIt.GoToBegin();
    b0ImageIt.GoToBegin();

    Eigen::VectorXd B(m_NumImages);
    Eigen::VectorXd D(7);
    CoefficientMatrixType weightMatrix = m_BMatrix;
    double resiMatrix;
    typename NumericTraits< ReferencePixelType >::AccumulateType b0 =
             NumericTraits< ReferencePixelType >::Zero;    
    // if a mask is present, iterate through mask image and skip zero voxels
    typename MaskSpatialObjectType::Pointer maskSpatialObject;
    if(this->m_MaskImagePresent)
    {
        maskSpatialObject =
                static_cast<MaskSpatialObjectType *>(this->ProcessObject::GetInput(1));
    }
    bool useMask(maskSpatialObject.IsNotNull());
    // Two cases here .
    // 1. If the Gradients have been specified in multiple images, we will create
    // 'n' iterators for each of the gradient images and solve the Stejskal-Tanner
    // equations for every pixel.
    // 2. If the Gradients have been specified in a single multi-component image,
    // one iterator will suffice to do the same.
    if ( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
       typedef ImageRegionConstIterator< GradientImageType > GradientIteratorType;
        std::vector< GradientIteratorType * > gradientItContainer;
                typename GradientImageType::Pointer gradientImagePointer =
                dynamic_cast< GradientImageType * >( this->ProcessObject::GetInput(0) );
            if(gradientImagePointer.IsNull())
            {
                itkExceptionMacro(<< "Invalid dynamic_cast");
            }                
            GradientIteratorType *git = new GradientIteratorType(
                        gradientImagePointer, outputRegionForThread);
            git->GoToBegin();
            gradientItContainer.push_back(git);
        for ( unsigned int i = 1; i <= m_NumImages; i++ )
        {
            gradientImagePointer =
                    dynamic_cast< GradientImageType * >( this->ProcessObject::GetInput(i+1) );
            if(gradientImagePointer.IsNull())
            {
                itkExceptionMacro(<< "Invalid dynamic_cast");
            }
            GradientIteratorType *git = new GradientIteratorType(
                        gradientImagePointer, outputRegionForThread);
            git->GoToBegin();
            gradientItContainer.push_back(git);
        }
        std::vector< double > pixelValue;    
        // Iterate over the reference and gradient images and solve the steskal
        // equations to reconstruct the Diffusion tensor.
        // See splweb.bwh.harvard.edu:8000/pages/papers/westin/ISMRM2002.pdf
        // "A Dual Tensor Basis Solution to the Stejskal-Tanner Equations for
        // DT-MRI"
        while ( !git->IsAtEnd() )
        {
            TensorPixelType tensor(0.0);
            //
            // if a mask is present, and we don't have a zero pixel
            // look up the voxel in the mask image corresponding to
            // the location of the current index.
            bool unmaskedPixel(true);
            if(useMask)
            {
                typename ImageRegionConstIteratorWithIndex<ReferenceImageType>::IndexType
                        index = git->GetIndex();
                typename ReferenceImageType::PointType point;
                gradientImagePointer->TransformIndexToPhysicalPoint(index,point);
                unmaskedPixel = maskSpatialObject->IsInside(point);
            }
            if ( unmaskedPixel )
            {
                for ( unsigned int i = 0; i < m_NumImages; i++ )
                {
                    GradientPixelType b = gradientItContainer[i]->Get();
                    if ( b == 0 )
                    {
                        B[i] = 0;
                        pixelValue.push_back(0);
                    }
                    else
                    {
                        B[i] = vcl_log( static_cast< double >( b ) );// LLS method
                        pixelValue.push_back(b);                        
                        if(this->m_Method == WLLS){
                            for (int j=0; j<7 ;j++) // Ln(S0) should not be weighted
                                weightMatrix(j,i) = m_BMatrix(j,i) * vcl_log( static_cast< double >(b) ) *
                                        vcl_log( static_cast< double >(b) ); //WLLS method
                        }
                    }
                    ++( *gradientItContainer[i] );
                }
                if ( m_NumberOfGradientDirections > 6 )
                {
                    m_TensorBasis = weightMatrix * m_BMatrix.transpose();
                }
                else
                {
                    m_TensorBasis = m_BMatrix.transpose();
                }
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_TensorBasis, Eigen::ComputeThinU | Eigen::ComputeThinV);
                if ( m_NumberOfGradientDirections > 6 )
                {
                    D = svd.solve(weightMatrix * B);
                }
                else
                {
                    D = svd.solve(B);
                }
		b0 = vcl_exp( D[0] );
                tensor(0, 0) = D[1];
                tensor(0, 1) = D[2];
                tensor(0, 2) = D[3];
                tensor(1, 1) = D[4];
                tensor(1, 2) = D[5];
                tensor(2, 2) = D[6];
                if(m_ADCImage){
		    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                    	if(this->m_BValue[i] != 0 && b0 != 0)
				B[i] = (vcl_log( static_cast< double >( b0 ) ) - B[i]) / this->m_BValue[i];
			else
				B[i] = 0;
		    }
                    ADCImageIt.Set(B.sum());
		}
                if(m_residImage){
                    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                        resiMatrix = m_BMatrix(1, i) * D[1] + m_BMatrix(2, i) * D[2] + m_BMatrix(3, i) * D[3] + 
				     m_BMatrix(4, i) * D[4] + m_BMatrix(5, i) * D[5] + m_BMatrix(6, i) * D[6];
                        B[i] = std::pow(static_cast< double >( pixelValue[i] ) - static_cast< double >( b0 ) *
                                vcl_exp( resiMatrix ), 2);
                    }
                    resiImageIt.Set(B.sum()/m_NumImages);
                }
                pixelValue.clear();                                
            }
            else
            {
                for ( unsigned int i = 0; i < m_NumberOfGradientDirections; i++ )
                {
                    ++( *gradientItContainer[i] );
                }
            }
            oit.Set(tensor);
            ++oit; // Output (reconstructed tensor image) iterator            
	    b0ImageIt.Set( b0 );
	    ++b0ImageIt;	    
            ++resiImageIt; // Output image of mean DWI residuals
            ++ADCImageIt; // ADC image output
        }
        for ( unsigned int i = 0; i < gradientItContainer.size(); i++ )
        {
            delete gradientItContainer[i];
        }
    }
    // The gradients are specified in a single multi-component image
    else if ( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
        typedef ImageRegionConstIterator< GradientImagesType >
                GradientIteratorType;
        typedef typename GradientImagesType::PixelType
                GradientVectorType;
        typename GradientImagesType::Pointer gradientImagePointer = NULL;
        // Would have liked a dynamic_cast here, but seems SGI doesn't like it
        // The enum will ensure that an inappropriate cast is not done
        gradientImagePointer = itkDynamicCastInDebugMode< GradientImagesType * >
                (this->ProcessObject::GetInput(0) );
        GradientIteratorType git(gradientImagePointer, outputRegionForThread);
        git.GoToBegin();

        while ( !git.IsAtEnd() )
        {
            GradientVectorType b = git.Get();
            TensorPixelType tensor(0.0);
            // if a mask is present, and we don't have a zero pixel
            // look up the voxel in the mask image corresponding to
            // the location of the current index.
            bool unmaskedPixel(true);
            if(useMask)
            {
                typename ImageRegionConstIteratorWithIndex<ReferenceImageType>::IndexType
                        index = git.GetIndex();
                typename ReferenceImageType::PointType point;
                gradientImagePointer->TransformIndexToPhysicalPoint(index,point);
                unmaskedPixel = maskSpatialObject->IsInside(point);
            }
            if ( unmaskedPixel )
            {
                for ( unsigned int i = 0; i < m_NumImages; i++ )
                {
                    if ( b[i] == 0 )
                    {
                        B[i] = 0;
                    }
                    else
                    {
                        B[i] = vcl_log( static_cast< double >( b[i] ) );// LLS method
                        if(this->m_Method == WLLS){
                            for (int j=0; j<7 ;j++) // Ln(S0) should not be weighted
                                weightMatrix(j, i) = m_BMatrix(j, i) * vcl_log( static_cast< double >(b[i]) ) *
                                        vcl_log( static_cast< double >(b[i]) ); //WLLS method
                        }
                    }
                }

                if ( m_NumberOfGradientDirections > 6 )
                {
                    m_TensorBasis = weightMatrix * m_BMatrix.transpose();
                }
                else
                {
                    m_TensorBasis = m_BMatrix.transpose();
                }

                //vnl_svd< double > pseudoInverseSolver(m_TensorBasis);
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(m_TensorBasis, Eigen::ComputeThinU | Eigen::ComputeThinV);
                if ( m_NumberOfGradientDirections > 6 )
                {
                    D = svd.solve(weightMatrix * B);
                }
                else
                {
                    D = svd.solve(B);
                }
				b0 = vcl_exp( D[0] );
                tensor(0, 0) = D[1];
                tensor(0, 1) = D[2];
                tensor(0, 2) = D[3];
                tensor(1, 1) = D[4];
                tensor(1, 2) = D[5];
                tensor(2, 2) = D[6];

                if(m_ADCImage){
		    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                    	if(this->m_BValue[i] != 0 && b0 != 0)
				B[i] = (vcl_log( static_cast< double >( b0 ) ) - B[i]) / this->m_BValue[i];
			else
				B[i] = 0;
		    }
                    ADCImageIt.Set(B.sum());
		}
                if(m_residImage){
                    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                        resiMatrix = m_BMatrix(1, i) * D[1] + m_BMatrix(2, i) * D[2] + m_BMatrix(3, i) * D[3] + 
				     m_BMatrix(4, i) * D[4] + m_BMatrix(5, i) * D[5] + m_BMatrix(6, i) * D[6];
                        B[i] = std::pow(static_cast< double >( b[i] ) - static_cast< double >( b0 ) *
                                vcl_exp( resiMatrix ), 2);
                    }
                    resiImageIt.Set(B.sum()/m_NumImages);
                }
            }
            oit.Set(tensor);
	    b0ImageIt.Set( b0 );
            ++oit; // Output (reconstructed tensor image) iterator
            ++git; // Gradient image iterator
            ++resiImageIt; // Output image of mean DWI residuals
            ++ADCImageIt; // ADC image output
	    ++b0ImageIt;
        }
    }
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::AfterThreadedGenerateData(){

}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::ComputeTensorBasis()
{
    if ( m_NumberOfGradientDirections < 6 )
    {
        itkExceptionMacro(<< "Not enough gradient directions supplied. Need to supply at least 6");
    }
    // This is only important if we are using a vector image. For
    // images added one at a time, this is not needed but doesn't hurt.
    std::vector< unsigned int > diffimageind;
    for ( typename GradientDirectionContainerType::ConstIterator gdcit = this->m_GradientDirectionContainer->Begin();
          gdcit != this->m_GradientDirectionContainer->End(); ++gdcit )
    {
            diffimageind.push_back( gdcit.Index() );
    }
    m_BMatrix.resize(m_NumImages, 7);
    for ( unsigned int m = 0; m < m_NumImages; m++ )
    {
		m_BMatrix(m,0) = 1; // Used for Ln(S0) estimate
        m_BMatrix(m,1) = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0] * this->m_BValue[m];
        m_BMatrix(m,2) = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1] * this->m_BValue[m];
        m_BMatrix(m,3) = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2] * this->m_BValue[m];
        m_BMatrix(m,4) = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1] * this->m_BValue[m];
        m_BMatrix(m,5) = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2] * this->m_BValue[m];
        m_BMatrix(m,6) = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2] * this->m_BValue[m];
    }
    if ( m_NumberOfGradientDirections > 6 )
    {
        m_TensorBasis = m_BMatrix.transpose() * m_BMatrix;
    }
    else
    {
        m_TensorBasis = m_BMatrix;
    }
    m_BMatrix.transposeInPlace();
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
const Image<TGradientImagePixelType,3> *
TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::GetGradientImage(unsigned index) const
{
    if ( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
        itkExceptionMacro(<< "Cannot retrieve individual gradient Image if "
                          << "all gradients are in a single image.");
    }
    // input 0 is either the single gradient image, or the reference
    // image. input 1 is either null or a mask image.
    return itkDynamicCastInDebugMode< const GradientImageType * >
            (this->ProcessObject::GetInput(index+2));
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::AddGradientImage(const GradientDirectionType & gradientDirection,
                   const GradientImageType *gradientImage)
{
    // Make sure crazy users did not call both AddGradientImage and
    // SetGradientImage
    if ( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
        itkExceptionMacro(<< "Cannot call both methods:"
                          << "AddGradientImage and SetGradientImage. Please call only one of them.");
    }
    // If the container to hold the gradient directions hasn't been allocated
    // yet, allocate it.
    if ( !this->m_GradientDirectionContainer )
    {
        this->m_GradientDirectionContainer = GradientDirectionContainerType::New();
    }
    m_GradientDirectionContainer->InsertElement(
                m_NumberOfGradientDirections, gradientDirection / gradientDirection.two_norm() );
    ++m_NumberOfGradientDirections;
    this->ProcessObject::SetNthInput( m_NumberOfGradientDirections+1,
                                      const_cast< GradientImageType * >( gradientImage ) );
    m_GradientImageTypeEnumeration = GradientIsInManyImages;
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::SetGradientImage(GradientDirectionContainerType *gradientDirection,
                   const GradientImagesType *gradientImage)
{
    // Make sure crazy users did not call both AddGradientImage and
    // SetGradientImage
    if ( m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        itkExceptionMacro(<< "Cannot call both methods:"
                          << "AddGradientImage and SetGradientImage. Please call only one of them.");
    }
    this->m_GradientDirectionContainer = gradientDirection;
    m_NumImages = gradientDirection->Size();

    this->m_NumberOfBaselineImages = 0;
    int i=0;
    newBValue.clear();
    if(m_BValue.size() == 0)
        std::cout<< "BValue should be defined before SetGradientImage()"<<std::endl;
    for ( typename GradientDirectionContainerType::Iterator it = this->m_GradientDirectionContainer->Begin();
          it != this->m_GradientDirectionContainer->End(); it++ )
    {
        if ( it.Value().one_norm() <= 0.0 )
        {
            this->m_NumberOfBaselineImages++;
            i++;
        }
        else // Normalize non-zero gradient directions
        {
            it.Value() = it.Value() / it.Value().two_norm();
            newBValue.push_back(m_BValue[i]); // fill non zero bvalues in a new vector
            i++;
        }
    }
    this->m_NumberOfGradientDirections = m_NumImages - this->m_NumberOfBaselineImages;
    //std::cout<<"number of baseline "<<m_NumberOfBaselineImages<<" others "<<m_NumberOfGradientDirections<<std::endl;
    // ensure that the gradient image we received has as many components as
    // the number of gradient directions
    if ( gradientImage->GetVectorLength() != this->m_NumberOfBaselineImages + this->m_NumberOfGradientDirections )
    {
        itkExceptionMacro(<< this->m_NumberOfGradientDirections << " gradients + " << this->m_NumberOfBaselineImages
                          << "baselines = " << this->m_NumberOfGradientDirections + this->m_NumberOfBaselineImages
                          << " directions specified but image has " << gradientImage->GetVectorLength()
                          << " components.");
    }
//    this->SetInput(gradientImage);
    this->ProcessObject::SetNthInput( 0,
                                      const_cast< GradientImagesType * >( gradientImage ) );
    m_GradientImageTypeEnumeration = GradientIsInASingleImage;
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void
TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::SetMaskSpatialObject(MaskSpatialObjectType *maskSpatialObject)
{
    this->ProcessObject::SetNthInput(1,maskSpatialObject);
    this->m_MaskImagePresent = true;
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void
TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::SetMaskImage(MaskImageType *maskImage)
{
    typename ImageMaskSpatialObject<3>::Pointer maskSpatialObject =
            ImageMaskSpatialObject<3>::New();
    maskSpatialObject->SetImage(maskImage);
    this->SetMaskSpatialObject(maskSpatialObject.GetPointer());
}

template< class TReferenceImagePixelType,
          class TGradientImagePixelType, class TTensorPixelType,
          class TMaskImageType >
void TensorReconstructionFilter_tbbnd< TReferenceImagePixelType,
TGradientImagePixelType,
TTensorPixelType,
TMaskImageType >
::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << "TensorBasisMatrix: " << m_TensorBasis << std::endl;
    os << indent << "Coeffs: " << m_BMatrix << std::endl;
    if ( m_GradientDirectionContainer )
    {
        os << indent << "GradientDirectionContainer: "
           << m_GradientDirectionContainer << std::endl;
    }
    else
    {
        os << indent
           << "GradientDirectionContainer: (Gradient directions not set)" << std::endl;
    }
    os << indent << "NumberOfGradientDirections: "
       << m_NumberOfGradientDirections << std::endl;
    os << indent << "NumberOfBaselineImages: "
       << m_NumberOfBaselineImages << std::endl;
    os << indent << "Threshold for reference B0 image: " << m_Threshold << std::endl;
    //os << indent << "BValue: " << m_BValue << std::endl;
    if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        os << indent << "Gradient images haven been supplied " << std::endl;
    }
    else if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        os << indent << "A multicomponent gradient image has been supplied" << std::endl;
    }
}
}  // end of namespace itk

#endif // __TENSORRECONSTRUCTIONFILTER_TBBND_TXX
