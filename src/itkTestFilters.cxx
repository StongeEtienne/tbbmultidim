
#include "itk/itkRicianNoiseCorrectionFilter_itk.txx"
#include "itk/itkTensorReconstructionFilter_itk.txx"
#include "itk/itkPSDTensorEstimationFilter_itk.txx"

#include "itknd/itkRicianNoiseCorrectionFilter_itknd.txx"
#include "itknd/itkTensorReconstructionFilter_itknd.txx"
#include "itknd/itkPSDTensorEstimationFilter_itknd.txx"

#include "tbb/itkRicianNoiseCorrectionFilter_tbb.txx"
#include "tbb/itkTensorReconstructionFilter_tbb.txx"
#include "tbb/itkPSDTensorEstimationFilter_tbb.txx"

#include "tbbnd/itkRicianNoiseCorrectionFilter_tbbnd.txx"
#include "tbbnd/itkTensorReconstructionFilter_tbbnd.txx"
#include "tbbnd/itkPSDTensorEstimationFilter_tbbnd.txx"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "crlDWIStudy.h"
#include "PositiveTensorFunctor.h"
#include <tclap/CmdLine.h>
#include "itkImage.h"
#include "tbb/tick_count.h"
/************************************************************************************************
 * \class	PSDTensorEstimationFilterTest
 *
 * \brief	PSDTensorEstimationFilterTest
 *
 * \author	Amir Jaberzadeh
 * \date	July 2015
*************************************************************************************************/

int main( int argc, char *argv[] )
{

    std::string *fixedImageFile = new std::string("");
    float thresh = 0;
    unsigned short thread = 0;
    unsigned short Maxiter = 400;
    short imageFilter = 0;
    std::string *outputImageFile = new std::string("");
    std::string *outputImage2File = new std::string("");
    std::string *initialTensor = new std::string("");
    std::string *resiImage = new std::string("");
    std::string *b0Image = new std::string("");
    std::string *methodName = new std::string("");
    double tolerance = 1.e-4;
    float variance = 0;
    try
    {
        TCLAP::CmdLine cmd("Computational Radiology Laboratory",' ' );
        TCLAP::UnlabeledValueArg<std::string> inputImageArg("fixedImage","Fixed Image File Name",true,"","InputImage",cmd);
        TCLAP::ValueArg<float> threshArg("t","threshold","threshold of DWI intensities",false,thresh,"Threshold",cmd);
        TCLAP::ValueArg<unsigned int> ThreadArg("c","NCores","No of cores used",false,thread,"NCores",cmd);
        TCLAP::ValueArg<unsigned int> MaxiterArg("r","MaxIteration","Maximum number of iterations",false,Maxiter,"MaxIteration",cmd);
        TCLAP::ValueArg<unsigned int> ImageFilterArg("i","ImageFilter","0:all, 1:ITK, 2:ITKND, 3:TBB, 4:TBBND",false,imageFilter,"ImageFilter",cmd);      
        TCLAP::UnlabeledValueArg<std::string> outputImageArg("outputImage","Output tensor image",true,"","OutputTensorImage.nrrd",cmd);
        TCLAP::UnlabeledValueArg<std::string> outputImage2Arg("numReps","Image showing no of iterations in each voxel",false,"","OptimizerIterationsImage.nrrd",cmd);
        TCLAP::ValueArg<std::string> methodNameArg("m","methodName","Name of method",false,"","Method name(CLLS, CWLLS, CNLS Or CWNLS)",cmd);
        TCLAP::ValueArg<std::string> TensorInput("p","iniTensor",	"Initial tensor used for residual Image",	false, "", "InputTensorImage.nrrd", cmd);
        TCLAP::ValueArg<std::string> resiOutput("o","resiImage",	"residual Image filename",	false, "", "OutputResidImage.nrrd", cmd);
        TCLAP::ValueArg<std::string> b0Output("b","b0Image",	"b0 Image filename",	false, "", "InputB0Image.nrrd", cmd);      
        TCLAP::ValueArg<float> VarianceArg("v","variance","Background noise variance",false,variance,"Noise variance",cmd); 
        cmd.parse(argc,argv);

        if (inputImageArg.isSet()) fixedImageFile = new std::string(inputImageArg.getValue());
        if (threshArg.isSet()) thresh = threshArg.getValue();
        if (ThreadArg.isSet()) thread = ThreadArg.getValue();
        if (MaxiterArg.isSet()) Maxiter = MaxiterArg.getValue();
        if (ImageFilterArg.isSet()) imageFilter = ImageFilterArg.getValue();
        if (outputImageArg.isSet()) outputImageFile = new std::string(outputImageArg.getValue());
        if (outputImage2Arg.isSet()) outputImage2File = new std::string(outputImage2Arg.getValue());
        if (TensorInput.isSet()) initialTensor = new std::string(TensorInput.getValue());
        if (resiOutput.isSet()) resiImage = new std::string(resiOutput.getValue());
        if (b0Output.isSet()) b0Image = new std::string(b0Output.getValue());      
        if (methodNameArg.isSet()) methodName = new std::string(methodNameArg.getValue());
        if (VarianceArg.isSet()) variance = VarianceArg.getValue();   
    }
    catch (TCLAP::ArgException& e)
    {
      std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
      exit(1);
    }
    typedef double PixelType;
    typedef itk::DiffusionTensor3D< double > TensorPixelType;
    typedef itk::VectorImage<PixelType, 3> VectorImageType;
    typedef VectorImageType::Pointer ImagePointerType;
    typedef itk::Image<PixelType, 3> ImageType;
    typedef crl::DWIStudy<PixelType> StudyType;
    StudyType inputStudy;
    inputStudy.LoadStudy(*fixedImageFile);

    tbb::tick_count starttime, endtime; // for timing
    double tbb_time;
    ImagePointerType img = inputStudy.DWIData;

    // RUN FOR ITK
    if(imageFilter == 1 || imageFilter == 0){

        // RicianNoiseCorrectionFilter
        typedef itk::RicianNoiseCorrectionFilter_itk< VectorImageType, ImageType, VectorImageType> RicianNoiseCorrectionFilterType;
        RicianNoiseCorrectionFilterType::Pointer ricianFilter = RicianNoiseCorrectionFilterType::New();
        ricianFilter->SetInput(img);
        ricianFilter->SetVariance(variance);
        ricianFilter->SetCoordinateTolerance(0.0001);
        ricianFilter->SetDirectionTolerance(0.0001);

        if(thread != 0){
            ricianFilter->SetNumberOfThreads( thread  );
        	//ricianFilter->SetManualThreads();
        }
        starttime = tbb::tick_count::now();
        ricianFilter->Update();
        endtime = tbb::tick_count::now();
        tbb_time = (endtime - starttime).seconds();

        std::cout << " ITK RicianNoiseCorrectionFilter time for "<<thread<<" threads = " << tbb_time -  ricianFilter->m_beforetime<< " seconds" << std::endl;

        // init tensor type
        typedef itk::TensorReconstructionFilter_itk< PixelType, PixelType, PixelType > TensorReconstructionFilterType;
        TensorReconstructionFilterType::Pointer tensorReconstructionFilter = TensorReconstructionFilterType::New();
        TensorReconstructionFilterType::GradientDirectionContainerType::Pointer
        BVecs = TensorReconstructionFilterType::GradientDirectionContainerType::New();
        for(int i = 0; i < inputStudy.Gradients.size(); ++i){
            BVecs->InsertElement(i, inputStudy.Gradients.at(i));
        }

        // TensorReconstructionFilter
		tensorReconstructionFilter->SetBValue(inputStudy.BValues);
        tensorReconstructionFilter->SetGradientImage(BVecs, img );
        tensorReconstructionFilter->SetMethod(TensorReconstructionFilterType::LLS);
        if(thread != 0){
        	tensorReconstructionFilter->SetNumberOfThreads( thread  );
        	//tensorReconstructionFilter->SetManualThreads();
        }
        tensorReconstructionFilter->SetDirectionTolerance(tolerance);
        tensorReconstructionFilter->SetCoordinateTolerance(tolerance);  
      
      	starttime = tbb::tick_count::now();
        tensorReconstructionFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();
      	std::cout << " ITK TensorReconstructionFilter time for " << thread << " threads = " 
                  << tbb_time - tensorReconstructionFilter->m_beforetime<< " seconds" << std::endl;

        // Correcting non positive semi definite tensors
        typedef TensorReconstructionFilterType::OutputImageType TensorType;        
        typedef itk::PositiveTensorImageFilter<TensorType, TensorType > CorrectionFilterType;
        CorrectionFilterType::Pointer correctionFilter = CorrectionFilterType::New();
        
        if(initialTensor->empty()){
            correctionFilter->SetInput(tensorReconstructionFilter->GetOutput());
            correctionFilter->Update();
        }
        typedef itk::ImageFileReader< TensorType > TensorReaderType;
        TensorReaderType::Pointer reader = TensorReaderType::New();
        if(!initialTensor->empty()){
            reader->SetFileName( *initialTensor );
            reader->Update();
        }
        // -------------------------------------------------------------------------
        // Run PSDTensorEstimation class
        typedef itk::PSDTensorEstimationFilter_itk< PixelType, ImageType, PixelType > PSDTensorEstimationFilterType;
        PSDTensorEstimationFilterType::Pointer PSDtensorEstimationFilter = PSDTensorEstimationFilterType::New();
        PSDtensorEstimationFilter->SetBValue(inputStudy.BValues);
        PSDtensorEstimationFilter->SetGradientImage( BVecs, img );
        PSDtensorEstimationFilter->SetDirectionTolerance(tolerance);
        PSDtensorEstimationFilter->SetCoordinateTolerance(tolerance);
           
	    if (methodName->empty())
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CLLS);
	    else if (methodName->compare("CWLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWLLS);
	    else if (methodName->compare("CWNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWNLS);

        if(thread != 0){
        	PSDtensorEstimationFilter->SetNumberOfThreads(thread);
        	//PSDtensorEstimationFilter->SetManualThreads();
        }
        if(!resiImage->empty()){
            PSDtensorEstimationFilter->SetMeanResiduals();
        }

        PSDtensorEstimationFilter->SetMaxIterations(Maxiter);
        PSDtensorEstimationFilter->SetThreshold( thresh );
        if(!initialTensor->empty()){
            PSDtensorEstimationFilter->SetInitialTensorImage(reader->GetOutput());
        }
        else{
            PSDtensorEstimationFilter->SetInitialTensorImage(correctionFilter->GetOutput());
        }
                
        PSDtensorEstimationFilter->SetInitialB0Image(tensorReconstructionFilter->GetOutput4());
        PSDtensorEstimationFilter->SetMeanResiduals();

      	starttime = tbb::tick_count::now();
        PSDtensorEstimationFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();

      	std::cout << " ITK PSDTensorEstimationFilter time for "<<thread<<" threads = " 
                  << tbb_time -  PSDtensorEstimationFilter->m_beforetime<< " seconds" << std::endl;
        // -------------------------------------------------------------------------
        // Write out the image of iterations. This code snippet goes to show that you
        // can use itk::ImageFileWriter to write an image of iterations.
        //
        typedef itk::ImageFileWriter<PSDTensorEstimationFilterType::OutputImageType1 > TensorWriterType1;
        TensorWriterType1::Pointer iterationWriter1 = TensorWriterType1::New();
        iterationWriter1->SetFileName( *outputImageFile );
        iterationWriter1->SetInput( PSDtensorEstimationFilter->GetOutput1() );
        iterationWriter1->Update();

	    if(!outputImage2File->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType2 > TensorWriterType2;
            TensorWriterType2::Pointer iterationWriter2 = TensorWriterType2::New();
            iterationWriter2->SetFileName( *outputImage2File );
            iterationWriter2->SetInput( PSDtensorEstimationFilter->GetOutput2() );
            iterationWriter2->Update();
	    }
        if(!resiImage->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType3 > TensorWriterType3;
            TensorWriterType3::Pointer iterationWriter3 = TensorWriterType3::New();
            iterationWriter3->SetFileName( *resiImage );
            iterationWriter3->SetInput( PSDtensorEstimationFilter->GetOutput3() );
            iterationWriter3->Update();
        }
        if(!b0Image->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType4 > TensorWriterType4;
            TensorWriterType4::Pointer iterationWriter4 = TensorWriterType4::New();
            iterationWriter4->SetFileName( *b0Image );
            iterationWriter4->SetInput( PSDtensorEstimationFilter->GetOutput4() );
            iterationWriter4->Update();
        }
    }


    // RUN FOR ITKND
    if(imageFilter == 2 || imageFilter == 0){

        // RicianNoiseCorrectionFilter
        typedef itk::RicianNoiseCorrectionFilter_itknd< VectorImageType, ImageType, VectorImageType> RicianNoiseCorrectionFilterType;
        RicianNoiseCorrectionFilterType::Pointer ricianFilter = RicianNoiseCorrectionFilterType::New();
        ricianFilter->SetInput(img);
        ricianFilter->SetVariance(variance);
        ricianFilter->SetCoordinateTolerance(0.0001);
        ricianFilter->SetDirectionTolerance(0.0001);

        if(thread != 0){
            ricianFilter->SetNumberOfThreads( thread  );
        	//ricianFilter->SetManualThreads();
        }
        starttime = tbb::tick_count::now();
        ricianFilter->Update();
        endtime = tbb::tick_count::now();
        tbb_time = (endtime - starttime).seconds();

        std::cout << " ITKND RicianNoiseCorrectionFilter time for "<<thread<<" threads = " << tbb_time -  ricianFilter->m_beforetime<< " seconds" << std::endl;

        // init tensor type
        typedef itk::TensorReconstructionFilter_itknd< PixelType, PixelType, PixelType > TensorReconstructionFilterType;
        TensorReconstructionFilterType::Pointer tensorReconstructionFilter = TensorReconstructionFilterType::New();
        TensorReconstructionFilterType::GradientDirectionContainerType::Pointer
        BVecs = TensorReconstructionFilterType::GradientDirectionContainerType::New();
        for(int i = 0; i < inputStudy.Gradients.size(); ++i){
            BVecs->InsertElement(i, inputStudy.Gradients.at(i));
        }

        // TensorReconstructionFilter
		tensorReconstructionFilter->SetBValue(inputStudy.BValues);
        tensorReconstructionFilter->SetGradientImage(BVecs, img );
        tensorReconstructionFilter->SetMethod(TensorReconstructionFilterType::LLS);
        if(thread != 0){
        	tensorReconstructionFilter->SetNumberOfThreads( thread  );
        	//tensorReconstructionFilter->SetManualThreads();
        }
        tensorReconstructionFilter->SetDirectionTolerance(tolerance);
        tensorReconstructionFilter->SetCoordinateTolerance(tolerance);  
      
      	starttime = tbb::tick_count::now();
        tensorReconstructionFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();
      	std::cout << " ITKND TensorReconstructionFilter time for " << thread << " threads = " 
                  << tbb_time - tensorReconstructionFilter->m_beforetime<< " seconds" << std::endl;

        // Correcting non positive semi definite tensors
        typedef TensorReconstructionFilterType::OutputImageType TensorType;        
        typedef itk::PositiveTensorImageFilter<TensorType, TensorType > CorrectionFilterType;
        CorrectionFilterType::Pointer correctionFilter = CorrectionFilterType::New();
        
        if(initialTensor->empty()){
            correctionFilter->SetInput(tensorReconstructionFilter->GetOutput());
            correctionFilter->Update();
        }
        typedef itk::ImageFileReader< TensorType > TensorReaderType;
        TensorReaderType::Pointer reader = TensorReaderType::New();
        if(!initialTensor->empty()){
            reader->SetFileName( *initialTensor );
            reader->Update();
        }
        // -------------------------------------------------------------------------
        // Run PSDTensorEstimation class
        typedef itk::PSDTensorEstimationFilter_itknd<PixelType, ImageType, PixelType > PSDTensorEstimationFilterType;
        PSDTensorEstimationFilterType::Pointer PSDtensorEstimationFilter = PSDTensorEstimationFilterType::New();
        PSDtensorEstimationFilter->SetBValue(inputStudy.BValues);
        PSDtensorEstimationFilter->SetGradientImage( BVecs, img );
        PSDtensorEstimationFilter->SetDirectionTolerance(tolerance);
        PSDtensorEstimationFilter->SetCoordinateTolerance(tolerance);
           
	    if (methodName->empty())
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CLLS);
	    else if (methodName->compare("CWLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWLLS);
	    else if (methodName->compare("CWNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWNLS);

        if(thread != 0){
        	PSDtensorEstimationFilter->SetNumberOfThreads(thread);
        	//PSDtensorEstimationFilter->SetManualThreads();
        }
        if(!resiImage->empty()){
            PSDtensorEstimationFilter->SetMeanResiduals();
        }

        PSDtensorEstimationFilter->SetMaxIterations(Maxiter);
        PSDtensorEstimationFilter->SetThreshold( thresh );
        if(!initialTensor->empty()){
            PSDtensorEstimationFilter->SetInitialTensorImage(reader->GetOutput());
        }
        else{
            PSDtensorEstimationFilter->SetInitialTensorImage(correctionFilter->GetOutput());
        }
                
        PSDtensorEstimationFilter->SetInitialB0Image(tensorReconstructionFilter->GetOutput4());
        PSDtensorEstimationFilter->SetMeanResiduals();

      	starttime = tbb::tick_count::now();
        PSDtensorEstimationFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();

      	std::cout << " ITKND PSDtensorEstimationFilter time for "<<thread<<" threads = " 
                  << tbb_time -  PSDtensorEstimationFilter->m_beforetime<< " seconds" << std::endl;
        // -------------------------------------------------------------------------
        // Write out the image of iterations. This code snippet goes to show that you
        // can use itk::ImageFileWriter to write an image of iterations.
        //
        typedef itk::ImageFileWriter<PSDTensorEstimationFilterType::OutputImageType1 > TensorWriterType1;
        TensorWriterType1::Pointer iterationWriter1 = TensorWriterType1::New();
        iterationWriter1->SetFileName( *outputImageFile );
        iterationWriter1->SetInput( PSDtensorEstimationFilter->GetOutput1() );
        iterationWriter1->Update();

	    if(!outputImage2File->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType2 > TensorWriterType2;
            TensorWriterType2::Pointer iterationWriter2 = TensorWriterType2::New();
            iterationWriter2->SetFileName( *outputImage2File );
            iterationWriter2->SetInput( PSDtensorEstimationFilter->GetOutput2() );
            iterationWriter2->Update();
	    }
        if(!resiImage->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType3 > TensorWriterType3;
            TensorWriterType3::Pointer iterationWriter3 = TensorWriterType3::New();
            iterationWriter3->SetFileName( *resiImage );
            iterationWriter3->SetInput( PSDtensorEstimationFilter->GetOutput3() );
            iterationWriter3->Update();
        }
        if(!b0Image->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType4 > TensorWriterType4;
            TensorWriterType4::Pointer iterationWriter4 = TensorWriterType4::New();
            iterationWriter4->SetFileName( *b0Image );
            iterationWriter4->SetInput( PSDtensorEstimationFilter->GetOutput4() );
            iterationWriter4->Update();
        }
    }

    // RUN FOR TBB
    if(imageFilter == 3 || imageFilter == 0){

        // RicianNoiseCorrectionFilter
        typedef itk::RicianNoiseCorrectionFilter_tbb< VectorImageType, ImageType, VectorImageType> RicianNoiseCorrectionFilterType;
        RicianNoiseCorrectionFilterType::Pointer ricianFilter = RicianNoiseCorrectionFilterType::New();
        ricianFilter->SetInput(img);
        ricianFilter->SetVariance(variance);
        ricianFilter->SetCoordinateTolerance(0.0001);
        ricianFilter->SetDirectionTolerance(0.0001);

        if(thread != 0){
            ricianFilter->SetNumberOfThreads( thread  );
        	ricianFilter->SetManualThreads();
        }
        starttime = tbb::tick_count::now();
        ricianFilter->Update();
        endtime = tbb::tick_count::now();
        tbb_time = (endtime - starttime).seconds();

        std::cout << " TBB RicianNoiseCorrectionFilter time for "<<thread<<" threads = " << tbb_time -  ricianFilter->m_beforetime<< " seconds" << std::endl;

        // init tensor type
        typedef itk::TensorReconstructionFilter_tbb<PixelType, PixelType, PixelType > TensorReconstructionFilterType;
        TensorReconstructionFilterType::Pointer tensorReconstructionFilter = TensorReconstructionFilterType::New();
        TensorReconstructionFilterType::GradientDirectionContainerType::Pointer
        BVecs = TensorReconstructionFilterType::GradientDirectionContainerType::New();
        for(int i = 0; i < inputStudy.Gradients.size(); ++i){
            BVecs->InsertElement(i, inputStudy.Gradients.at(i));
        }

        // init tensorReconstructionFilter
		tensorReconstructionFilter->SetBValue(inputStudy.BValues);
        tensorReconstructionFilter->SetGradientImage(BVecs, img );
        tensorReconstructionFilter->SetMethod(TensorReconstructionFilterType::LLS);
        if(thread != 0){
        	tensorReconstructionFilter->SetNumberOfThreads( thread  );
        	tensorReconstructionFilter->SetManualThreads();
        }
        tensorReconstructionFilter->SetDirectionTolerance(tolerance);
        tensorReconstructionFilter->SetCoordinateTolerance(tolerance);  
      

        //TIME
      	starttime = tbb::tick_count::now();
        tensorReconstructionFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();
      	std::cout << " TBB TensorReconstructionFilter time for " << thread << " threads = " 
                  << tbb_time - tensorReconstructionFilter->m_beforetime<< " seconds" << std::endl;

        // Correcting non positive semi definite tensors
        typedef TensorReconstructionFilterType::OutputImageType TensorType;        
        typedef itk::PositiveTensorImageFilter<TensorType, TensorType > CorrectionFilterType;
        CorrectionFilterType::Pointer correctionFilter = CorrectionFilterType::New();
        
        if(initialTensor->empty()){
            correctionFilter->SetInput(tensorReconstructionFilter->GetOutput());
            correctionFilter->Update();
        }
        typedef itk::ImageFileReader< TensorType > TensorReaderType;
        TensorReaderType::Pointer reader = TensorReaderType::New();
        if(!initialTensor->empty()){
            reader->SetFileName( *initialTensor );
            reader->Update();
        }
        // -------------------------------------------------------------------------
        // Run PSDTensorEstimation class
        typedef itk::PSDTensorEstimationFilter_tbb<PixelType, ImageType, PixelType > PSDTensorEstimationFilterType;
        PSDTensorEstimationFilterType::Pointer PSDtensorEstimationFilter = PSDTensorEstimationFilterType::New();
        PSDtensorEstimationFilter->SetBValue(inputStudy.BValues);
        PSDtensorEstimationFilter->SetGradientImage( BVecs, img );
        PSDtensorEstimationFilter->SetDirectionTolerance(tolerance);
        PSDtensorEstimationFilter->SetCoordinateTolerance(tolerance);
           
	    if (methodName->empty())
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CLLS);
	    else if (methodName->compare("CWLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWLLS);
	    else if (methodName->compare("CWNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWNLS);

        if(thread != 0){
        	PSDtensorEstimationFilter->SetNumberOfThreads(thread);
        	PSDtensorEstimationFilter->SetManualThreads();
        }
        if(!resiImage->empty()){
            PSDtensorEstimationFilter->SetMeanResiduals();
        }

        PSDtensorEstimationFilter->SetMaxIterations(Maxiter);
        PSDtensorEstimationFilter->SetThreshold( thresh );
        if(!initialTensor->empty()){
            PSDtensorEstimationFilter->SetInitialTensorImage(reader->GetOutput());
        }
        else{
            PSDtensorEstimationFilter->SetInitialTensorImage(correctionFilter->GetOutput());
        }
                
        PSDtensorEstimationFilter->SetInitialB0Image(tensorReconstructionFilter->GetOutput4());
        PSDtensorEstimationFilter->SetMeanResiduals();

      	starttime = tbb::tick_count::now();
        PSDtensorEstimationFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();

      	std::cout << " TBB PSDtensorEstimationFilter time for "<<thread<<" threads = " 
                  << tbb_time -  PSDtensorEstimationFilter->m_beforetime<< " seconds" << std::endl;
        // -------------------------------------------------------------------------
        // Write out the image of iterations. This code snippet goes to show that you
        // can use itk::ImageFileWriter to write an image of iterations.
        //
        typedef itk::ImageFileWriter<PSDTensorEstimationFilterType::OutputImageType1 > TensorWriterType1;
        TensorWriterType1::Pointer iterationWriter1 = TensorWriterType1::New();
        iterationWriter1->SetFileName( *outputImageFile );
        iterationWriter1->SetInput( PSDtensorEstimationFilter->GetOutput1() );
        iterationWriter1->Update();

	    if(!outputImage2File->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType2 > TensorWriterType2;
            TensorWriterType2::Pointer iterationWriter2 = TensorWriterType2::New();
            iterationWriter2->SetFileName( *outputImage2File );
            iterationWriter2->SetInput( PSDtensorEstimationFilter->GetOutput2() );
            iterationWriter2->Update();
	    }
        if(!resiImage->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType3 > TensorWriterType3;
            TensorWriterType3::Pointer iterationWriter3 = TensorWriterType3::New();
            iterationWriter3->SetFileName( *resiImage );
            iterationWriter3->SetInput( PSDtensorEstimationFilter->GetOutput3() );
            iterationWriter3->Update();
        }
        if(!b0Image->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType4 > TensorWriterType4;
            TensorWriterType4::Pointer iterationWriter4 = TensorWriterType4::New();
            iterationWriter4->SetFileName( *b0Image );
            iterationWriter4->SetInput( PSDtensorEstimationFilter->GetOutput4() );
            iterationWriter4->Update();
        }
    }
    
    // RUN FOR TBBND
    if(imageFilter == 4 || imageFilter == 0){

        // RicianNoiseCorrectionFilter
        typedef itk::RicianNoiseCorrectionFilter_tbbnd< VectorImageType, ImageType, VectorImageType> RicianNoiseCorrectionFilterType;
        RicianNoiseCorrectionFilterType::Pointer ricianFilter = RicianNoiseCorrectionFilterType::New();
        ricianFilter->SetInput(img);
        ricianFilter->SetVariance(variance);
        ricianFilter->SetCoordinateTolerance(0.0001);
        ricianFilter->SetDirectionTolerance(0.0001);

        if(thread != 0){
            ricianFilter->SetNumberOfThreads( thread  );
        	//ricianFilter->SetManualThreads();
        }
        starttime = tbb::tick_count::now();
        ricianFilter->Update();
        endtime = tbb::tick_count::now();
        tbb_time = (endtime - starttime).seconds();

        std::cout << " TBBND RicianNoiseCorrectionFilter time for "<<thread<<" threads = " << tbb_time -  ricianFilter->m_beforetime<< " seconds" << std::endl;

        // init tensor type
        typedef itk::TensorReconstructionFilter_tbbnd< PixelType, PixelType, PixelType > TensorReconstructionFilterType;
        TensorReconstructionFilterType::Pointer tensorReconstructionFilter = TensorReconstructionFilterType::New();
        TensorReconstructionFilterType::GradientDirectionContainerType::Pointer
        BVecs = TensorReconstructionFilterType::GradientDirectionContainerType::New();
        for(int i = 0; i < inputStudy.Gradients.size(); ++i){
            BVecs->InsertElement(i, inputStudy.Gradients.at(i));
        }

        // TensorReconstructionFilter
		tensorReconstructionFilter->SetBValue(inputStudy.BValues);
        tensorReconstructionFilter->SetGradientImage(BVecs, img );
        tensorReconstructionFilter->SetMethod(TensorReconstructionFilterType::LLS);
        if(thread != 0){
        	tensorReconstructionFilter->SetNumberOfThreads( thread  );
        	//tensorReconstructionFilter->SetManualThreads();
        }
        tensorReconstructionFilter->SetDirectionTolerance(tolerance);
        tensorReconstructionFilter->SetCoordinateTolerance(tolerance);  
      
      	starttime = tbb::tick_count::now();
        tensorReconstructionFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();
      	std::cout << " TBBND TensorReconstructionFilter time for " << thread << " threads = " 
                  << tbb_time - tensorReconstructionFilter->m_beforetime<< " seconds" << std::endl;

        // Correcting non positive semi definite tensors
        typedef TensorReconstructionFilterType::OutputImageType TensorType;        
        typedef itk::PositiveTensorImageFilter<TensorType, TensorType > CorrectionFilterType;
        CorrectionFilterType::Pointer correctionFilter = CorrectionFilterType::New();
        
        if(initialTensor->empty()){
            correctionFilter->SetInput(tensorReconstructionFilter->GetOutput());
            correctionFilter->Update();
        }
        typedef itk::ImageFileReader< TensorType > TensorReaderType;
        TensorReaderType::Pointer reader = TensorReaderType::New();
        if(!initialTensor->empty()){
            reader->SetFileName( *initialTensor );
            reader->Update();
        }
        // -------------------------------------------------------------------------
        // Run PSDTensorEstimation class
        typedef itk::PSDTensorEstimationFilter_tbbnd<PixelType, ImageType, PixelType > PSDTensorEstimationFilterType;
        PSDTensorEstimationFilterType::Pointer PSDtensorEstimationFilter = PSDTensorEstimationFilterType::New();
        PSDtensorEstimationFilter->SetBValue(inputStudy.BValues);
        PSDtensorEstimationFilter->SetGradientImage( BVecs, img );
        PSDtensorEstimationFilter->SetDirectionTolerance(tolerance);
        PSDtensorEstimationFilter->SetCoordinateTolerance(tolerance);
           
	    if (methodName->empty())
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CNLS);
	    else if (methodName->compare("CLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CLLS);
	    else if (methodName->compare("CWLLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWLLS);
	    else if (methodName->compare("CWNLS") == 0)
            	PSDtensorEstimationFilter->SetMethod(PSDTensorEstimationFilterType::CWNLS);

        if(thread != 0){
        	PSDtensorEstimationFilter->SetNumberOfThreads(thread);
        	//PSDtensorEstimationFilter->SetManualThreads();
        }
        if(!resiImage->empty()){
            PSDtensorEstimationFilter->SetMeanResiduals();
        }

        PSDtensorEstimationFilter->SetMaxIterations(Maxiter);
        PSDtensorEstimationFilter->SetThreshold( thresh );
        if(!initialTensor->empty()){
            PSDtensorEstimationFilter->SetInitialTensorImage(reader->GetOutput());
        }
        else{
            PSDtensorEstimationFilter->SetInitialTensorImage(correctionFilter->GetOutput());
        }
                
        PSDtensorEstimationFilter->SetInitialB0Image(tensorReconstructionFilter->GetOutput4());
        PSDtensorEstimationFilter->SetMeanResiduals();

      	starttime = tbb::tick_count::now();
        PSDtensorEstimationFilter->Update();
      	endtime = tbb::tick_count::now();
      	tbb_time = (endtime - starttime).seconds();

      	std::cout << " TBBND PSDtensorEstimationFilter time for "<<thread<<" threads = " 
                  << tbb_time -  PSDtensorEstimationFilter->m_beforetime<< " seconds" << std::endl;
        // -------------------------------------------------------------------------
        // Write out the image of iterations. This code snippet goes to show that you
        // can use itk::ImageFileWriter to write an image of iterations.
        //
        typedef itk::ImageFileWriter<PSDTensorEstimationFilterType::OutputImageType1 > TensorWriterType1;
        TensorWriterType1::Pointer iterationWriter1 = TensorWriterType1::New();
        iterationWriter1->SetFileName( *outputImageFile );
        iterationWriter1->SetInput( PSDtensorEstimationFilter->GetOutput1() );
        iterationWriter1->Update();

	    if(!outputImage2File->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType2 > TensorWriterType2;
            TensorWriterType2::Pointer iterationWriter2 = TensorWriterType2::New();
            iterationWriter2->SetFileName( *outputImage2File );
            iterationWriter2->SetInput( PSDtensorEstimationFilter->GetOutput2() );
            iterationWriter2->Update();
	    }
        if(!resiImage->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType3 > TensorWriterType3;
            TensorWriterType3::Pointer iterationWriter3 = TensorWriterType3::New();
            iterationWriter3->SetFileName( *resiImage );
            iterationWriter3->SetInput( PSDtensorEstimationFilter->GetOutput3() );
            iterationWriter3->Update();
        }
        if(!b0Image->empty()){
            typedef itk::ImageFileWriter< PSDTensorEstimationFilterType::OutputImageType4 > TensorWriterType4;
            TensorWriterType4::Pointer iterationWriter4 = TensorWriterType4::New();
            iterationWriter4->SetFileName( *b0Image );
            iterationWriter4->SetInput( PSDtensorEstimationFilter->GetOutput4() );
            iterationWriter4->Update();
        }


    }

}

