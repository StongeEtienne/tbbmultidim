#ifndef __ITKPSDTENSORESTIMATIONFILTER_ITKND_TXX
#define __ITKPSDTENSORESTIMATIONFILTER_ITKND_TXX

 
#include "itkCRLImageToImageFilter.h"
#include "itkCRLImageToImageFilter.txx"

#include "itkPSDTensorEstimationFilter_itknd.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkArray.h"
#include "itkImageMaskSpatialObject.h"
#include "vnl/vnl_vector.h"

#define USE_ITK_OPTIMIZER
#ifdef USE_ITK_OPTIMIZER
#include "itkPowellOptimizer.h"
#include <itkAmoebaOptimizer.h>
#include <itkParticleSwarmOptimizer.h>
#else
#include "itkNLOPTOptimizers.h"
#endif

#include <cfloat>
#include "tbb/tick_count.h"


/************************************************************************************************
 * \class	PSDTensorEstimationFilter_itknd
 *
 * \brief	PSDTensorEstimationFilter_itknd
 *
 * \author	Amir Jaberzadeh, Etienne St-Onge
 * \date	June 2017
*************************************************************************************************/

namespace itk
{
using namespace std;
template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
TTensorPixelType,
TMaskImageType >
::PSDTensorEstimationFilter_itknd()
{
    // At least 1 inputs is necessary for a vector image.
    // For images added one at a time we need at least six
    this->SetNumberOfRequiredInputs(0);
    this->SetNumberOfRequiredOutputs(4);
    this->SetNthOutput( 0, ( OutputImageType1::New() ).GetPointer());
    this->SetNthOutput( 1, ( OutputImageType2::New() ).GetPointer());
    this->SetNthOutput( 2, ( OutputImageType3::New() ).GetPointer());
    this->SetNthOutput( 3, ( OutputImageType4::New() ).GetPointer());
    m_NumberOfGradientDirections = 0;
    m_NumberOfBaselineImages = 1;
    m_NumImages = 1;
    m_Threshold = NumericTraits< float >::min();
    m_GradientImageTypeEnumeration = Else;
    m_GradientDirectionContainer = NULL;
    m_TensorBasis.set_identity();
    m_InitialTensor = NULL;
    m_B0Image = NULL;
    m_MaskImagePresent = false;
    m_TensorImagePresent = false;
    m_B0ImagePresent = false;
    m_Method = CNLS;
    m_Maxiter = 500;
    m_MeanResiduals = false;
    m_AllResiduals = false;
}


template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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
//    residImage = OutputImageType2::New();
//    residImage->SetSpacing(this->GetOutput2()->GetSpacing());
//    residImage->SetOrigin(this->GetOutput2()->GetOrigin());
//    residImage->SetDirection(this->GetOutput2()->GetDirection());
//    residImage->SetRegions(this->GetOutput2()->GetLargestPossibleRegion());
//    residImage->Allocate();
//    residImage->FillBuffer(0);
    // If there's a mask, make sure it matches the dimensions of the
    // gradient image(s).

    endtime = tbb::tick_count::now();
    m_beforetime = (endtime - starttime).seconds();

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

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
TTensorPixelType,
TMaskImageType >
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
itk::ThreadIdType ThreadId)
{
	float twiceminnonzero = 2.*FLT_EPSILON;
    typename OutputImageType1::Pointer outputImage =
            this->GetOutput1();
    typename OutputImageType2::Pointer funcEvalImage =
            this->GetOutput2();
    typename OutputImageType3::Pointer residImage =
            this->GetOutput3();
    typename OutputImageType4::Pointer b0Image = 
    	    this->GetOutput4();        

    ImageRegionIterator< OutputImageType1 > oit(outputImage, outputRegionForThread);
    ImageRegionIterator< OutputImageType1 > tensorIterator(m_InitialTensor, outputRegionForThread);
    ImageRegionIterator< OutputImageType4 > B0Iterator(m_B0Image, outputRegionForThread);
    ImageRegionIterator< OutputImageType4 > OutputB0Iterator(b0Image, outputRegionForThread);    
    ImageRegionIterator< OutputImageType2 > funcEvalIterator(funcEvalImage, outputRegionForThread);
    ImageRegionIterator< OutputImageType3 > residImageIt(residImage, outputRegionForThread);

    residImageIt.GoToBegin();
    oit.GoToBegin();
    tensorIterator.GoToBegin();
    B0Iterator.GoToBegin();
    funcEvalIterator.GoToBegin();
    OutputB0Iterator.GoToBegin();
    
    vnl_vector< double > B(m_NumImages);
    vnl_vector< double > D(6);
    // if a mask is present, iterate through mask image and skip zero voxels
    typename MaskSpatialObjectType::Pointer maskSpatialObject;
    if(this->m_MaskImagePresent)
    {
        maskSpatialObject =
                static_cast<MaskSpatialObjectType *>(this->ProcessObject::GetInput(1));
    }
    bool useMask(maskSpatialObject.IsNotNull());
    
	#ifdef USE_ITK_OPTIMIZER
//         	typedef itk::AmoebaOptimizer OptimizerType; // ITK Optimizer
		typedef itk::PowellOptimizer OptimizerType;
//		typedef itk::ParticleSwarmOptimizer OptimizerType;		
	#else
        	typedef itk::NLOPTOptimizers OptimizerType;   // NLOPT Optimizer
	#endif

    // Create the cost function
    typedef itk::MyCostFunction2<TInputImagePixelType, OutputImageType3> CostFunType;
    typename CostFunType::Pointer costFun = CostFunType::New();
    costFun->SetMethod(this->m_Method);
    // Create the optimizer
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetCostFunction(costFun);


        // Set the parameters
	#ifdef USE_ITK_OPTIMIZER
		//itk::Optimizer::ScalesType scale(costFun->GetNumberOfParameters()); 
		//scale.Fill(1);
		//scale[6] = 1e3;
		//optimizer->SetScales(scale);
/*		
		optimizer->SetMaximalNumberOfIterations( m_Maxiter );     // Particle Swarm
		std::vector<double> lbound, ubound;
		for (int i=0; i<7; i++)
			{
			lbound.push_back(-9e-2);
			ubound.push_back(9e-2);
			}
    	ParticleSwarmOptimizer::ParameterBoundsType bounds;
    	for ( int i=0; i<lbound.size(); i++ )
      	{
      		std::pair<double,double> value;
      		value.first = lbound[i];
      		value.second = ubound[i];
      		bounds.push_back( value );
      	}
      	unsigned int numberOfParticles = 5;
		double xTolerance = 1e-3;
		double fTolerance = 1;		
		optimizer->SetParameterBounds( bounds );
 		optimizer->SetNumberOfParticles( numberOfParticles );
		optimizer->SetParametersConvergenceTolerance( xTolerance,
							costFun->GetNumberOfParameters() );
		optimizer->SetFunctionConvergenceTolerance( fTolerance );
*/

///*
        optimizer->SetMaximumIteration( m_Maxiter ); // Powell
 		optimizer->SetStepLength(1e-3);
		optimizer->SetStepTolerance(1e-4);
		optimizer->SetValueTolerance(1e-1);
//*/		
		/*optimizer->SetMaximumNumberOfIterations( m_Maxiter ); // AmoebaOptimizer
		optimizer->SetOptimizeWithRestarts(true);
		optimizer->SetFunctionConvergenceTolerance(1e-4); 
 		optimizer->SetParametersConvergenceTolerance(1e-4);*/
	#else
        	optimizer->SetMaxEval(m_Maxiter);
        	optimizer->SetAlgorithm(OptimizerType::NLOPT_LN_NEWUOA);
//        	optimizer->SetXTolRel(1e-3);
        	optimizer->SetXTolAbs(1e-4);
	#endif

    // Set the initial position
    OptimizerType::ParametersType  pinit(costFun->GetNumberOfParameters()),
    			     	   pBest(costFun->GetNumberOfParameters());        
    std::vector< double > pixelValue;
    double resiMatrix;
    double R0,R1,R2,R3,R4,R5;
    
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

        while ( !git->IsAtEnd() )
        {
            typename NumericTraits< ReferencePixelType >::AccumulateType b0 =
                    NumericTraits< ReferencePixelType >::Zero;
            if(m_B0ImagePresent)        
            	b0 = B0Iterator.Get();
            else
            	b0 = m_Threshold + 1;
            TensorPixelType tensor(0.0);
            tensor(0,0) = twiceminnonzero;
            tensor(1,1) = twiceminnonzero;
            tensor(2,2) = twiceminnonzero;
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
            if ( ( b0 != 0 ) && unmaskedPixel && ( b0 >= m_Threshold ) )
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
                        B[i] = vcl_log( static_cast< double >( b ));
                        pixelValue.push_back(static_cast< double >( b ));
                    }
                    ++( *gradientItContainer[i] );
                }
                try {
                    costFun->SetPixelLogValue(&B);
                    costFun->SetDesignMatrix(&m_BMatrix);
                    costFun->SetPixelIntensityValue(&pixelValue);
                    tensor = tensorIterator.Get();
//                    D[0] = tensor(0, 0);
//                    D[1] = tensor(0, 1);
//                    D[2] = tensor(0, 2);
//                    D[3] = tensor(1, 1);
//                    D[4] = tensor(1, 2);
//                    D[5] = tensor(2, 2);
/* ------------------ Tensors can be initialized with a given TensorImage.
*   These notations are taken from paper " Investigation of Anomalous Estimates of Tensor-Derived
*   Quantities in Diffusion Tensor Imaging" equations [4], [10]  ------------------------*/

                R0 = sqrt(tensor(0, 0) + 1e-8); // Dxx = R0^2 
                R5 = tensor(0, 2)/R0;// Dxz = R0R5;
                R3 = tensor(0, 1)/R0;// Dxy = R0R3  
                R1 = sqrt(tensor(1, 1) - R3*R3 + 1e-8);// Dyy = R1^2 + R3^2
                R4 = (tensor(1, 2) - R3*R5)/R1;// Dyz = R1R4 + R3R5
                R2 = sqrt(tensor(2, 2) - R4*R4 - R5*R5 + 1e-8);// Dzz = R2 ^2 + R4^2 + R5^2 
                            	
                            pinit[0] = R0;
                            pinit[1] = R1;
                            pinit[2] = R2;
                            pinit[3] = R3;
                            pinit[4] = R4;
                            pinit[5] = R5;
		    	    pinit[6] = b0 * 1e-5;   // S0 
		    	    
//   ----------- Switch next part with the previous one for initializing optimizer with known values
/*
                    pinit[0] = 1E-4;
                    pinit[1] = 1E-4;
                    pinit[2] = 1E-4;
                    pinit[3] = 0;
                    pinit[4] = 0;
                    pinit[5] = 0;
		    pinit[6] = b0 * 1e-5;  // S0
*/		    
// ------------------------------------
//		    optimizer->ClearSwarm(); // Particle Swarm
                    optimizer->SetInitialPosition(pinit);
                    optimizer->StartOptimization();
//                    std::cout<<"  "<< optimizer->GetCurrentCost();
//                    std::cout<< optimizer->GetNumberOfEvaluations() << " ";
//                    std::cout<<"Return code "<<(int)optimizer->GetErrorCode();
//                    std::cout<<" = "<<optimizer->GetErrorCodeDescription()<<endl;

		#ifdef USE_ITK_OPTIMIZER
//		    funcEvalIterator.Set(optimizer->GetCurrentIteration());
		#else
                    funcEvalIterator.Set(optimizer->GetNumberOfEvaluations());
                    if ( !optimizer->isSuccessful() )
                    {
                        std::cout<<"Fatal error"<<std::endl<<std::endl;
                    }
                    else
                    {		
		#endif
		
                        pBest = optimizer->GetCurrentPosition();
                        tensor(0, 0) = pBest[0]*pBest[0];
                        tensor(0, 1) = pBest[0]*pBest[3];
                        tensor(0, 2) = pBest[0]*pBest[5];
                        tensor(1, 1) = pBest[1]*pBest[1] + pBest[3]*pBest[3];
                        tensor(1, 2) = pBest[1]*pBest[4] + pBest[3]*pBest[5];
                        tensor(2, 2) = pBest[2]*pBest[2] + pBest[4]*pBest[4] + pBest[5]*pBest[5];
                        b0 = pBest[6] * 1e5;
 		#ifndef USE_ITK_OPTIMIZER
                    }
		#endif                                              
                    }
                /*-------------------------------------
                If catched an exception, show a message
                -------------------------------------*/
                catch (itk::ExceptionObject& e)
                {
                    cout<< "  Exception: "<<endl<< e.GetDescription();
                }
/*---------------- Uncomment this part to compute residuals based on Tensors calculated in this class using nonlinear
 * methods, otherwise the residuals are computed using user's given initial tensor image or the ones computed with LLS
 * method passed by PSDTensorEstimationFilterTest.cxx file. ------------------------ */
                D[0] = tensor(0, 0);
                D[1] = tensor(0, 1);
                D[2] = tensor(0, 2);
                D[3] = tensor(1, 1);
                D[4] = tensor(1, 2);
                D[5] = tensor(2, 2);
//  ---------------------------------- --------------------
                if(m_MeanResiduals){
                    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                        resiMatrix = m_BMatrix[0][i] * D[0] + m_BMatrix[1][i] * D[1] + m_BMatrix[2][i] * D[2] +
                                     m_BMatrix[3][i] * D[3] + m_BMatrix[4][i] * D[4] + m_BMatrix[5][i] * D[5];
                        B[i] = std::pow(static_cast< double >( pixelValue[i] ) - static_cast< double >( b0 ) * 
                        vcl_exp( resiMatrix ), 2);
                    }
                    residImageIt.Set(B.sum()/static_cast< double >(m_NumImages));                    
                }
                pixelValue.clear();
            }
            else
            {
                for ( unsigned int i = 0; i < m_NumImages; i++ )
                {
                    ++( *gradientItContainer[i] );
                }
            }
            ++tensorIterator;
            OutputB0Iterator.Set(b0);
            ++OutputB0Iterator;
            ++B0Iterator;
            ++funcEvalIterator;
            oit.Set(tensor);
            ++oit; // Output (reconstructed tensor image) iterator
            ++residImageIt; 
        }
        for ( unsigned int i = 0; i < gradientItContainer.size(); i++ )
        {
            delete gradientItContainer[i];
        }       
    }
    
    // The gradients are specified in a single multi-component image
    else if ( m_GradientImageTypeEnumeration == GradientIsInASingleImage )
    {
        typedef ImageRegionConstIteratorWithIndex< GradientImagesType >
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
            typename NumericTraits< ReferencePixelType >::AccumulateType b0 =
                    NumericTraits< ReferencePixelType >::Zero;
            if(m_B0ImagePresent)        
            	b0 = B0Iterator.Get();
            else
            	b0 = m_Threshold + 1;
            TensorPixelType tensor(0.0);
            tensor(0,0) = twiceminnonzero;
            tensor(1,1) = twiceminnonzero;
            tensor(2,2) = twiceminnonzero;
            //
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
            if ( ( b0 != 0 ) && unmaskedPixel && ( b0 >= m_Threshold ) )
            {
                for ( unsigned int i = 0; i < m_NumImages ; i++ )
                {
                    if ( b[i] == 0 )
                    {
                        B[i] = 0;
                        pixelValue.push_back(0);
                    }
                    else
                    {
                        B[i] = vcl_log( static_cast< double >( b[i] ));
                        pixelValue.push_back(static_cast< double >( b[i] ));
                    }
                }
                try {
                    costFun->SetPixelLogValue(&B);
                    costFun->SetDesignMatrix(&m_BMatrix);
                    costFun->SetPixelIntensityValue(&pixelValue);
                    tensor = tensorIterator.Get();
//                    D[0] = tensor(0, 0);
//                    D[1] = tensor(0, 1);
//                    D[2] = tensor(0, 2);
//                    D[3] = tensor(1, 1);
//                    D[4] = tensor(1, 2);
//                    D[5] = tensor(2, 2);
/* ------------------ Tensors can be initialized with a given TensorImage.
*   These notations are taken from paper " Investigation of Anomalous Estimates of Tensor-Derived
*   Quantities in Diffusion Tensor Imaging" equations [4], [10]  ------------------------*/

                R0 = sqrt(tensor(0, 0) + 1e-8); // Dxx = R0^2 
                R5 = tensor(0, 2)/R0;// Dxz = R0R5;
                R3 = tensor(0, 1)/R0;// Dxy = R0R3  
                R1 = sqrt(tensor(1, 1) - R3*R3 + 1e-8);// Dyy = R1^2 + R3^2
                R4 = (tensor(1, 2) - R3*R5)/R1;// Dyz = R1R4 + R3R5
                R2 = sqrt(tensor(2, 2) - R4*R4 - R5*R5 + 1e-8);// Dzz = R2 ^2 + R4^2 + R5^2 
                
                            pinit[0] = R0;
                            pinit[1] = R1;
                            pinit[2] = R2;
                            pinit[3] = R3;
                            pinit[4] = R4;
                            pinit[5] = R5;
		    	    pinit[6] = b0 * 1e-5;   // S0 
		    	    		    	    
//   ----------- Switch next part with the previous one for initializing optimizer with known values
/*
                    pinit[0] = 1E-2;
                    pinit[1] = 1E-2;
                    pinit[2] = 1E-2;
                    pinit[3] = 0;
                    pinit[4] = 0;
                    pinit[5] = 0;
		    pinit[6] = b0;  // S0
*/		      
// ------------------------------------
//		    optimizer->ClearSwarm(); // Particle Swarm
                    optimizer->SetInitialPosition(pinit);
                    optimizer->StartOptimization();
                    pixelValue.clear();

//                    std::cout<<"  "<< optimizer->GetCurrentCost();
//                    std::cout<< optimizer->GetNumberOfEvaluations() << " ";
//                    std::cout<<"Return code "<<(int)optimizer->GetErrorCode();
//                    std::cout<<" = "<<optimizer->GetErrorCodeDescription()<<endl;

		#ifdef USE_ITK_OPTIMIZER
//		    funcEvalIterator.Set(optimizer->GetCurrentIteration());
		#else
                    funcEvalIterator.Set(optimizer->GetNumberOfEvaluations());
                    if ( !optimizer->isSuccessful() )
                    {
                        std::cout<<"Fatal error"<<std::endl<<std::endl;
                    }
                    else
                    {		
		#endif
		
                        pBest = optimizer->GetCurrentPosition();
                        tensor(0, 0) = pBest[0]*pBest[0];
                        tensor(0, 1) = pBest[0]*pBest[3];
                        tensor(0, 2) = pBest[0]*pBest[5];
                        tensor(1, 1) = pBest[1]*pBest[1] + pBest[3]*pBest[3];
                        tensor(1, 2) = pBest[1]*pBest[4] + pBest[3]*pBest[5];
                        tensor(2, 2) = pBest[2]*pBest[2] + pBest[4]*pBest[4] + pBest[5]*pBest[5];
                        b0 = pBest[6] * 1e5;
 		#ifndef USE_ITK_OPTIMIZER
                    }
		#endif                         
                    }
                /*-------------------------------------
                If catched an exception, show a message
                -------------------------------------*/
                catch (itk::ExceptionObject& e)
                {
                    cout<< "  Exception: "<<endl<< e.GetDescription();
                }
/*---------------- Uncomment this part to compute residuals based on Tensors calculated in this class using nonlinear
 * methods, otherwise the residuals are computed using user's given initial tensor image or the ones computed with LLS
 * method passed by PSDTensorEstimationFilterTest.cxx file. ------------------------ */
                D[0] = tensor(0, 0);
                D[1] = tensor(0, 1);
                D[2] = tensor(0, 2);
                D[3] = tensor(1, 1);
                D[4] = tensor(1, 2);
                D[5] = tensor(2, 2);
//  ---------------------------------- --------------------
                if(m_MeanResiduals){
                    for ( unsigned int i = 0; i < m_NumImages; i++ )
                    {
                        resiMatrix = m_BMatrix[0][i] * D[0] + m_BMatrix[1][i] * D[1] + m_BMatrix[2][i] * D[2] +
                                     m_BMatrix[3][i] * D[3] + m_BMatrix[4][i] * D[4] + m_BMatrix[5][i] * D[5];
                        B[i] = std::pow(static_cast< double >( b[i] ) - static_cast< double >( b0 ) * 
                        vcl_exp( resiMatrix ), 2);
                    }
                    residImageIt.Set(B.sum()/static_cast< double >(m_NumImages));
                }
            }
	    
            ++tensorIterator;
            OutputB0Iterator.Set(b0);
            ++OutputB0Iterator;
            ++B0Iterator;
            ++funcEvalIterator;
            oit.Set(tensor);
            ++oit; // Output (reconstructed tensor image) iterator
            ++git; // Gradient image iterator
            ++residImageIt;
        }
    }
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
TTensorPixelType,
TMaskImageType >
::AfterThreadedGenerateData()
{
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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
    m_BMatrix.set_size(m_NumImages, 6);
    for ( unsigned int m = 0; m < m_NumImages; m++ )
    {
        m_BMatrix[m][0] = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0] * this->m_BValue[m];
        m_BMatrix[m][1] = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1] * this->m_BValue[m];
        m_BMatrix[m][2] = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[0]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2] * this->m_BValue[m];
        m_BMatrix[m][3] = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1] * this->m_BValue[m];
        m_BMatrix[m][4] = -2 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[1]
                * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2] * this->m_BValue[m];
        m_BMatrix[m][5] = -1 * m_GradientDirectionContainer->ElementAt(diffimageind[m])[2]
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
    m_BMatrix.inplace_transpose();
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
const Image<TInputImagePixelType,3> *
PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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
        this->m_GradientDirectionContainer = typename GradientDirectionContainerType::New();
    }
    m_GradientDirectionContainer->InsertElement(
                m_NumberOfGradientDirections, gradientDirection / gradientDirection.two_norm() );
    ++m_NumberOfGradientDirections;
    this->ProcessObject::SetNthInput( m_NumberOfGradientDirections+1,
                                      const_cast< GradientImageType * >( gradientImage ) );
    m_GradientImageTypeEnumeration = GradientIsInManyImages;
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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
    for ( typename GradientDirectionContainerType::Iterator it = this->m_GradientDirectionContainer->Begin();
          it != this->m_GradientDirectionContainer->End(); it++ )
    {
        if ( it.Value().one_norm() <= 0.0 )
        {
            this->m_NumberOfBaselineImages++;
        }
        else // Normalize non-zero gradient directions
        {
            it.Value() = it.Value() / it.Value().two_norm();
        }
    }
    this->m_NumberOfGradientDirections = m_NumImages - this->m_NumberOfBaselineImages;
    // ensure that the gradient image we received has as many components as
    // the number of gradient directions
    if ( gradientImage->GetVectorLength() != this->m_NumberOfBaselineImages + this->m_NumberOfGradientDirections )
    {
        itkExceptionMacro(<< this->m_NumberOfGradientDirections << " gradients + " << this->m_NumberOfBaselineImages
                          << "baselines = " << this->m_NumberOfGradientDirections + this->m_NumberOfBaselineImages
                          << " directions specified but image has " << gradientImage->GetVectorLength()
                          << " components.");
    }
    this->ProcessObject::SetNthInput( 0,
                                      const_cast< GradientImagesType * >( gradientImage ) );
    m_GradientImageTypeEnumeration = GradientIsInASingleImage;
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void
PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
TTensorPixelType,
TMaskImageType >
::SetMaskSpatialObject(MaskSpatialObjectType *maskSpatialObject)
{
    this->ProcessObject::SetNthInput(1,maskSpatialObject);
    this->m_MaskImagePresent = true;
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void
PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
TTensorPixelType,
TMaskImageType >
::SetMaskImage(MaskImageType *maskImage)
{
    typename ImageMaskSpatialObject<3>::Pointer maskSpatialObject =
            ImageMaskSpatialObject<3>::New();
    maskSpatialObject->SetImage(maskImage);
    this->SetMaskSpatialObject(maskSpatialObject.GetPointer());
}

template< class TInputImagePixelType,
          class TOutputImageType, class TTensorPixelType,
          class TMaskImageType >
void PSDTensorEstimationFilter_itknd< TInputImagePixelType,
TOutputImageType,
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
//    os << indent << "BValue: " << m_BValue << std::endl;
    if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        os << indent << "Gradient images haven been supplied " << std::endl;
    }
    else if ( this->m_GradientImageTypeEnumeration == GradientIsInManyImages )
    {
        os << indent << "A multicomponent gradient image has been supplied" << std::endl;
    }
}

//template< class TInputImagePixelType,
//          class TOutputImageType, class TTensorPixelType,
//          class TMaskImageType >
//DataObject::Pointer PSDTensorEstimationFilter_itknd< TInputImagePixelType,
//TOutputImageType,
//TTensorPixelType,
//TMaskImageType >
//::MakeOutput(unsigned int idx)
//{
//  DataObject::Pointer output;

//  switch ( idx )
//    {
//    case 0:
//      output = ( OutputImageType1::New() ).GetPointer();
//      break;
//    case 1:
//      output = ( OutputImageType2::New() ).GetPointer();
//      break;
//    default:
//      std::cerr << "No output " << idx << std::endl;
//      output = NULL;
//      break;
//    }
//  return output.GetPointer();
//}

} // END Namespace itk
#endif // __ITKPSDTENSORESTIMATIONFILTER_ITK_TXX
