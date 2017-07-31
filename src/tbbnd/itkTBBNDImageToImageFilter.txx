#ifndef __itkTBBNDImageToImageFilter_hxx
#define __itkTBBNDImageToImageFilter_hxx

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    // We need that for TBB to correctly detect the Windows version
    //#define WINVER _WIN32_WINNT_WIN7
    //#define _WIN32_WINNT _WIN32_WINNT_WIN7
#include "windows.h"
#endif

#include "itkTBBNDImageToImageFilter.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/compat/thread"

#include "itkImageSource.h"
#include "itkImageRegionSplitterBase.h"
#include "itkOutputDataObjectIterator.h"
#include "vnl/vnl_math.h"

#define JOB_PER_THREAD_RATIO 15

/************************************************************************************************
 * \class    TBBNDImageToImageFilter
 *
 * \brief    TBBNDImageToImageFilter.
 *
 * \author    Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
 * \date    September 2016
*************************************************************************************************/

namespace itk {


    /**********************************************************************************************//**
     * \class    TBBNDFunctor
     *
     * \brief    TBB functor to execute jobs in parallel.
     *
     * \author    Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
     *
     * \tparam    TInputImage     Type of the input image.
     * \tparam    TOutputImage    Type of the output image.
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    class TBBNDFunctor 
    {
    public:
        typedef TBBNDFunctor        Self;
        typedef TOutputImage    OutputImageType;
        typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
        typedef typename TOutputImage::SizeType OutputImageSizeType;
        typedef typename OutputImageType::RegionType OutputImageRegionType;

        itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
        itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

        typedef TBBNDImageToImageFilter<TInputImage,TOutputImage> TbbImageFilterType;

        TBBNDFunctor(TbbImageFilterType *tbbFilter, const OutputImageSizeType& outputSize):
        m_TBBFilter(tbbFilter), m_OutputSize(outputSize) {}

        void operator() ( const tbb::blocked_range<int>& r ) const
        {
            typename TOutputImage::SizeType size = m_OutputSize;
            typename TOutputImage::IndexType index;
            index.Fill(0);

            if (m_TBBFilter->GetNbReduceDimensions() > 0)
            {
                unsigned int i = OutputImageDimension - (unsigned int)m_TBBFilter->GetNbReduceDimensions();

                index[i] = r.begin();
                size[i] = 1;
                while (i < OutputImageDimension - 1)
                {
                    index[i+1] = index[i] / m_OutputSize[i];
                    index[i] = index[i] % m_OutputSize[i];
                    size[i+1] = 1;
                    i++;
                }
            }

            // Construct an itk::ImageRegion
            OutputImageRegionType myRegion(index, size);

            // Run the TBBGenerateData method! (equivalent of ThreadedGenerateData)
            m_TBBFilter->TBBGenerateData(myRegion);
        }

    private:
        TbbImageFilterType *m_TBBFilter;
        OutputImageSizeType m_OutputSize;
    };


    // Constructor
    template< typename TInputImage, typename TOutputImage >
    TBBNDImageToImageFilter< TInputImage, TOutputImage >::TBBNDImageToImageFilter()
    {
        // By default, do not define the number of threads.
        // Let TBB doing that.
        this->SetNumberOfThreads(0);

        // By default, Automatic NbReduceDimensions
        this->SetNbReduceDimensions(-1);
    }

    // Destructor
    template< typename TInputImage, typename TOutputImage >
    TBBNDImageToImageFilter< TInputImage, TOutputImage >::~TBBNDImageToImageFilter()
    {}

    /************************************************************************************************
     * \fn    template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage > ::GenerateData()
     *
     * \brief    New default implementation for GenerateData() to use TBB
     *
     * \author  Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
     *
     * \tparam    typename TInputImage     Type of the typename t input image.
     * \tparam    typename TOutputImage    Type of the typename t output image.
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::GenerateData()
    {
        // Get the size of the requested region
        typename TOutputImage::ConstPointer output = static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
        typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();

        // Call a method that can be overriden by a subclass to allocate
        // memory for the filter's outputs
        this->AllocateOutputs();

        // Call a method that can be overridden by a subclass to perform
        // some calculations prior to splitting the main computations into
        // separate threads
        this->BeforeThreadedGenerateData();

        // Set up the number of threads with default
        // if it was not previously set
        if (this->GetNumberOfThreads() <= 0)
        {
            this->SetNumberOfThreads(tbb::task_scheduler_init::default_num_threads());
        }

        // Set up the TBB task_scheduler
        tbb::task_scheduler_init tbb_init(this->GetNumberOfThreads());

        // Generate the number of Jobs
        // based on the OutputImageDimension, NumberOfThreads and NbReduceDimensions
        this->GenerateNumberOfJobs();

        // Debug Output
        //std::cout << "TBB: " << this->GetNumberOfJobs() << "jobs, " << this->GetNumberOfThreads() << "threads;" << std::endl;

        // Do the task decomposition using parallel_for
        tbb::parallel_for( 
            tbb::blocked_range<int>(0, this->GetNumberOfJobs()),
            TBBNDFunctor<TInputImage, TOutputImage>(this, outputSize),
            tbb::simple_partitioner()); // force size 1 step

        // Call a method that can be overridden by a subclass to perform
        // some calculations after all the threads have completed
        this->AfterThreadedGenerateData();
    }


    template< typename TInputImage, typename TOutputImage >
    const ThreadIdType & TBBNDImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfThreads() const
    {
        return m_TBBNumberOfThreads;
    }

    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::SetNumberOfThreads(ThreadIdType nbThreads)
    {
        this->m_TBBNumberOfThreads = nbThreads;
    }


    /************************************************************************************************
     * \fn    template< typename TInputImage, typename TOutputImage > unsigned int TBBNDImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
     *
     * \brief    Gets the number of jobs.
     *             
     *             \warning This function only returns a valid value when called from AllocateOutputs,
     *             BeforeThreadedGenerateData(), TBBGenerateData() and AfterThreadedGenerateData(),
     *             ie when the input image is known
     *
     * \author  Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
     *
     * \tparam    TInputImage     Type of the input image.
     * \tparam    TOutputImage    Type of the output image.
     *
     * \return    The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    unsigned int TBBNDImageToImageFilter< TInputImage, TOutputImage >::GetNumberOfJobs() const
    {
        return m_TBBNumberOfJobs;
    }


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage >::SetNumberOfJobs(unsigned int nbJobs)
     *
     * \brief   Sets the number of jobs (Internal).
     *
     * \author  Amir Jaberzadeh, Benoit Scherrer and Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::SetNumberOfJobs(unsigned int nbJobs)
    {
        this->m_TBBNumberOfJobs = nbJobs;
    }



    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > unsigned int TBBNDImageToImageFilter< TInputImage, TOutputImage >::GetNbReduceDimensions() const
     *
     * \brief   Gets the number of dimension to separate for the Jobs multithreading
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    unsigned int TBBNDImageToImageFilter< TInputImage, TOutputImage >::GetNbReduceDimensions() const
    {
        return m_TBBNbReduceDimensions;
    }


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage >::SetNbReduceDimensions(int nbReduceDim)
     *
     * \brief   Set the number of dimension to separate and multithread each section.
     *          (nbReduceDim < 0  : negative number for automatic splitting)
     *
     *          \example : for a 3D image (volume) with the shape 30x10x5
     *          nbReduceDim == 0  : Will generate a single (1) Job with the whole image (size 30x10x5)
     *          nbReduceDim == 1  : Will generate 5 Jobs with the slices (size 30x10)
     *          nbReduceDim == 2  : Will generate 50 Jobs with the lines (size 30)
     *          nbReduceDim == 3  : Will generate 1500 Jobs with each voxel (size 1)
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::SetNbReduceDimensions(int nbReduceDim)
    {
        if (nbReduceDim > (int)OutputImageDimension)
        {
            this->m_TBBNbReduceDimensions = (int)OutputImageDimension;
        }
        else
        {
            this->m_TBBNbReduceDimensions = nbReduceDim;
        }
    }


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage >::GenerateNumberOfJobs()
     *
     * \brief   Generate the number Jobs based on the NbReduceDimensions
     *              or based on the NumberOfThreads and the Image Dimension (if NbReduceDimensions was not set).
     *
     *          \warning  This function must be called after the NumberOfThreads is set.
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::GenerateNumberOfJobs()
    {
        // Get the size of the requested region
        typename TOutputImage::ConstPointer output = 
        	static_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
        typename TOutputImage::SizeType outputSize = output->GetRequestedRegion().GetSize();
        const int imgDim = OutputImageDimension;
        const int nbReduceDim = GetNbReduceDimensions();

        // Generate the number of job
        if (nbReduceDim < 0)
        {
            // Heuristic to automatically determines m_TBBNbReduceDimensions
            m_TBBNbReduceDimensions = 0;
            m_TBBNumberOfJobs = 1;
            int currentDim = imgDim - 1;

            // Minimum Number of Jobs, based on the Number of thread
            unsigned int minNbJobs = JOB_PER_THREAD_RATIO * GetNumberOfThreads();
            while( currentDim >= 0 && m_TBBNumberOfJobs < minNbJobs )
            {
                ++m_TBBNbReduceDimensions;
                m_TBBNumberOfJobs *= outputSize[currentDim];
                --currentDim;
            }
        }
        else
        {
            // If manually chosen m_TBBNbReduceDimensions
            assert(nbReduceDim <= imgDim);

            m_TBBNumberOfJobs = 1;
            for (int i = imgDim - nbReduceDim; i < imgDim; ++i)
            {
                m_TBBNumberOfJobs *= outputSize[i];
            }
        }
    }



    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
     *
     * \brief   Use *TBBGenerateData()* instead of ThreadedGenerateData with TBBNDImageToImageFilter
     *
     *          \warning TBBNDImageToImageFilter doesn't support threadId
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
    {
        std::ostringstream message;
        message << "itk::ERROR: " << this->GetNameOfClass() << "(" << this << "): "
            << "Use *TBBGenerateData()* instead of ThreadedGenerateData() with TBBNDImageToImageFilter \n"
            << " TBBNDImageToImageFilter doesn't  support ThreadId ";

	std::cout << message.str() << std::endl;
        ExceptionObject e_(__FILE__, __LINE__, message.str().c_str(),ITK_LOCATION);
        throw e_;
    }


    /************************************************************************************************
     * \fn  template< typename TInputImage, typename TOutputImage > void TBBNDImageToImageFilter< TInputImage, TOutputImage >::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId)
     *
     * \brief   If an imaging filter can be implemented as a TBB multithreaded algorithm,
     *          the filter will provide an implementation of TBBGenerateData().
     *          This superclass will automatically split the output image into a number of pieces,
     *          spawn multiple threads, and call TBBGenerateData() in each thread.
     *          Prior to spawning threads, the BeforeThreadedGenerateData() method is called.
     *          After all the threads have completed, the AfterThreadedGenerateData() method is called.
     *
     *          \warning TBBNDImageToImageFilter doesn't support threadId
     *
     * \author  Etienne St-Onge
     *
     * \tparam  TInputImage     Type of the input image.
     * \tparam  TOutputImage    Type of the output image.
     *
     * \return  The number of jobs
     **************************************************************************************************/
    template< typename TInputImage, typename TOutputImage >
    void TBBNDImageToImageFilter< TInputImage, TOutputImage >::TBBGenerateData(const OutputImageRegionType& outputRegionForThread)
    {
        std::ostringstream message;
        message << "itk::ERROR: " << this->GetNameOfClass()
            << "(" << this << "): " << "Subclass should override this method!!!";
        ExceptionObject e_(__FILE__, __LINE__, message.str().c_str(),ITK_LOCATION);
        throw e_;
    }


}  //namespace itk

#endif
