/**********************************************************************************************//**
 * \file	crlMFMIO.cxx
 *
 * \brief	The crl::MFMIO class. 
*************************************************************************************************/

/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#include "crlMFMIO.h"

#include "configuration.h"
#include "crlFileName.h"
#include "crlMFMUtils.h"

#include <itkImageFileReader.h>
#include "itkImageFileWriter.h"

#include <itkDOMNodeXMLWriter.h>
#include <itkDOMNodeXMLReader.h>
#include <itkFancyString.h>

using namespace std;

namespace crl {		

	/**********************************************************************************************//**
	 * \fn	MFMIO::MFMIO()
	 *
	 * \brief	Constructor.
	 *
	 * \author	Benoit Scherrer
	 * \date	August 2013
	 **************************************************************************************************/
	MFMIO::MFMIO()
	{
	
	}

	/**********************************************************************************************//**
	 * \fn	MFMIO::~MFMIO()
	 *
	 * \brief	Destructor.
	 *
	 * \author	Benoit Scherrer
	 * \date	August 2013
	 **************************************************************************************************/
	MFMIO::~MFMIO()
	{

	}

	/**********************************************************************************************//**
	 * \fn	void MFMIO::ReadFromCommandLine( const std::vector<std::string>& fileNames,
	 * 		const std::string& fractionName );
	 *
	 * \brief	Reads from command line. The preferred way is to read a single XML file 
	 * 			containing the links to all the image files (tensor, fractions, B=0, etc...). 
	 * 			This is what is done when fileNames.size()==1 and the file extension is xml or mfm.
	 * 			However, for compatibility reasons, this function can also read the tensors and 
	 * 			fractions from the filenames in fileNames and fractionName.
	 *
	 * \author	Benoit Scherrer
	 * \date	August 2013
	 *
	 * \param	fileNames   	List of names of the files.
	 * \param	fractionName	Name of the fraction.
	 **************************************************************************************************/
	void MFMIO::ReadFromCommandLine( const std::vector<std::string>& fileNames, const std::string& fractionName, const std::string& kappaName )
	{
		bool readWithXML=false;

		if ( fileNames.size()==1 )
		{
			std::string fileext = crl::FileName(fileNames[0]).getExtension();
			std::transform(fileext.begin(), fileext.end(), fileext.begin(), ::tolower);

			if ( fileext=="mfm" || fileext=="xml" ) readWithXML = true;
		}

		if ( readWithXML )
		{
			std::cout<<"- Read DCI model (XML file)..."<<std::endl;
			ReadXML(fileNames[0]);
		}
		else
		{
			std::cout<<"- Read DCI model..."<<std::endl;

			MFM_LoadTensors<TensorImageType>(fileNames, m_TensorImages);
			if ( fractionName!="" )
			{
				std::cout<<"  " << fractionName.c_str() << std::endl;
				m_FractionsImage = MFM_OpenFractionsImageAsVectorImage<float>(fractionName);
			}
			else
			{
				// If no fraction, consider no isotropic and just equal fraction for all tensors
				m_FractionsImage = FractionsImageType::New();
				m_FractionsImage->SetRegions(m_TensorImages[0]->GetLargestPossibleRegion());
				m_FractionsImage->SetSpacing(m_TensorImages[0]->GetSpacing());
				m_FractionsImage->SetOrigin(m_TensorImages[0]->GetOrigin());
				m_FractionsImage->SetDirection(m_TensorImages[0]->GetDirection());
				m_FractionsImage->SetNumberOfComponentsPerPixel(m_TensorImages.size());
				m_FractionsImage->Allocate();

				FractionsImageType::PixelType initValue(m_TensorImages.size());
				initValue.Fill( 1.0 / ((double)m_TensorImages.size()));

				itk::ImageRegionIterator<FractionsImageType> it( m_FractionsImage, m_FractionsImage->GetLargestPossibleRegion());
				while ( !it.IsAtEnd() )
				{
					it.Set(initValue);
					++it;
				}
			}

			if ( kappaName!="" )
			{
				std::cout<<"  " << kappaName.c_str() << std::endl;
				m_GammaWidthImage = MFM_OpenFractionsImageAsVectorImage<float>(kappaName);		// use the same function to read the kappas
			}
			else
			{
				// If no fraction, consider no isotropic and just equal fraction for all tensors
				m_GammaWidthImage = KappasImageType::New();
				m_GammaWidthImage->SetRegions(m_TensorImages[0]->GetLargestPossibleRegion());
				m_GammaWidthImage->SetSpacing(m_TensorImages[0]->GetSpacing());
				m_GammaWidthImage->SetOrigin(m_TensorImages[0]->GetOrigin());
				m_GammaWidthImage->SetDirection(m_TensorImages[0]->GetDirection());
				m_GammaWidthImage->SetNumberOfComponentsPerPixel(m_TensorImages.size());
				m_GammaWidthImage->Allocate();

				KappasImageType::PixelType initValue(m_TensorImages.size());
				initValue.Fill( 999999999 );

				itk::ImageRegionIterator<FractionsImageType> it( m_GammaWidthImage, m_GammaWidthImage->GetLargestPossibleRegion());
				while ( !it.IsAtEnd() )
				{
					it.Set(initValue);
					++it;
				}
			}
		}

	}


	void MFMIO::clearTensorImages()
	{
		m_TensorImages.clear();
	}

	void MFMIO::addTensorImage(const TensorImageType::Pointer& image )
	{
		m_TensorImages.push_back(image);
	}

	void MFMIO::setFractionsImage(const FractionsImageType::Pointer& image )
	{
		m_FractionsImage = image;
	}

	void MFMIO::setKappasImage(const KappasImageType::Pointer& image )
	{
		m_GammaWidthImage = image;
	}

	void MFMIO::setB0Image(const FloatImageType::Pointer& image )
	{
		m_B0Image = image;
	}

	void MFMIO::setDisoImage(const FloatImageType::Pointer& image )
	{
		m_DisoImage = image;
	}

	const std::vector<MFMIO::TensorImageType::Pointer>& MFMIO::getTensorImages()  
	{
		return m_TensorImages ;
	}
	const MFMIO::TensorImageType::Pointer& MFMIO::getTensorImage(int id) 
	{
		return m_TensorImages.at(id);
	}


	const MFMIO::FractionsImageType::Pointer& MFMIO::getFractionsImage()  
	{
		return m_FractionsImage ;
	}
	const MFMIO::KappasImageType::Pointer&	MFMIO::getKappasImage()  
	{
		return m_GammaWidthImage ;
	}
	const MFMIO::FloatImageType::Pointer&	MFMIO::getB0Image()  
	{
		return m_B0Image ;
	}
	const MFMIO::FloatImageType::Pointer&	MFMIO::getDisoImage() 
	{
		return m_DisoImage ;
	}


	int MFMIO::getNumberOfTensors() const
	{
		return m_TensorImages.size();
	}
	
	bool MFMIO::freeDiffusion() const
	{
		if ( m_FractionsImage.IsNull() ) return false;
		else return  ( m_FractionsImage->GetNumberOfComponentsPerPixel() % getNumberOfTensors() == 1 );
	}

		
	void MFMIO::ReadXML ( const std::string& fileName )
	{
		std::cout<<"- Loading "<< fileName.c_str() <<std::endl;

		//-----------------------------------------
		// Read XML
		//-----------------------------------------
		itk::DOMNode::Pointer xmlRootNode;

		itk::DOMNodeXMLReader::Pointer reader = itk::DOMNodeXMLReader::New();
		reader->SetFileName( fileName );
		reader->Update();
		xmlRootNode = reader->GetOutput();

		//-----------------------------------------
		// Check
		//-----------------------------------------
		if ( xmlRootNode->GetName() != "MFM" )
			throw itk::ExceptionObject(__FILE__,__LINE__,"Invalid XML/MFM file.", "" );

		//-----------------------------------------
		// Load
		//-----------------------------------------
		loadTensors( xmlRootNode );
		loadFractions( xmlRootNode );
		loadB0( xmlRootNode );
		loadDiso( xmlRootNode );
		loadDIAMOND( xmlRootNode );		
	}

	void MFMIO::loadTensors( itk::DOMNode::Pointer node )
	{
		//-----------------------------------------
		// Get the XML node
		//-----------------------------------------
		itk::DOMNode* oNode = node->GetChild( "Tensors" );
		if ( oNode==NULL ) return;

		//-----------------------------------------
		// Read tensor files with ITK
		//-----------------------------------------
		std::cout<<"  - Loading the tensor file(s)..."<<std::endl;
		itk::DOMNode::ConstChildrenListType listTensors;
		oNode->GetAllChildren(listTensors);

		typedef itk::ImageFileReader<TensorImageType> TensorImageReaderType;
		TensorImageReaderType::Pointer reader = TensorImageReaderType::New();
		for ( unsigned int i=0; i<listTensors.size() ; i++ )
		{
			std::string fileName = listTensors.at(i)->GetAttribute("file" );
			std::cout<<"    " << fileName.c_str() <<std::endl;
			reader->SetFileName( fileName );
			reader->Update();
			m_TensorImages.push_back(reader->GetOutput());
			reader->GetOutput()->DisconnectPipeline();
		}
	}

	void MFMIO::loadFractions( itk::DOMNode::Pointer node )
	{
		//-----------------------------------------
		// Get the XML node
		//-----------------------------------------
		itk::DOMNode* oNode = node->GetChild( "Fractions" );
		if ( oNode==NULL ) return;

		//-----------------------------------------
		// Read file with ITK
		//-----------------------------------------
		std::string fileName = oNode->GetAttribute( "file" );
		std::cout<<"  - Loading the fractions from <" << fileName.c_str() << ">"<<std::endl;
		typedef itk::ImageFileReader<FractionsImageType> FractionsImageReaderType;
		FractionsImageReaderType::Pointer reader = FractionsImageReaderType::New();
		reader->SetFileName( fileName );
		reader->Update();
		
		this->m_FractionsImage = reader->GetOutput();
	}

	void MFMIO::loadB0( itk::DOMNode::Pointer node )
	{
		//-----------------------------------------
		// Get the XML node
		//-----------------------------------------
		itk::DOMNode* oNode = node->GetChild( "B0" );
		if ( oNode==NULL ) return;

		//-----------------------------------------
		// Read file with ITK
		//-----------------------------------------
		std::string fileName = oNode->GetAttribute( "file" );
		std::cout<<"  - Loading the B=0 image from <" << fileName.c_str() << ">"<<std::endl;
		typedef itk::ImageFileReader<FloatImageType> FloatImageReaderType;
		FloatImageReaderType::Pointer reader = FloatImageReaderType::New();
		reader->SetFileName( fileName );
		reader->Update();
		
		this->m_B0Image = reader->GetOutput();
	}

	void MFMIO::loadDiso( itk::DOMNode::Pointer node )
	{
		//-----------------------------------------
		// Get the XML node
		//-----------------------------------------
		itk::DOMNode* oNode = node->GetChild( "Diso" );
		if ( oNode==NULL ) return;

		//-----------------------------------------
		// Read file with ITK
		//-----------------------------------------
		std::string fileName = oNode->GetAttribute( "file" );
		std::cout<<"  - Loading the Diso image from <" << fileName.c_str() << ">"<<std::endl;
		typedef itk::ImageFileReader<FloatImageType> FloatImageReaderType;
		FloatImageReaderType::Pointer reader = FloatImageReaderType::New();
		reader->SetFileName( fileName );
		reader->Update();
		
		this->m_DisoImage = reader->GetOutput();
	}

	void MFMIO::loadDIAMOND( itk::DOMNode::Pointer node )
	{
		//-----------------------------------------
		// Get the XML node
		//-----------------------------------------
		itk::DOMNode* oNode = node->GetChild( "DIAMONDGamma" );
		if ( oNode==NULL ) return;

		//-----------------------------------------
		// Read file with ITK
		//-----------------------------------------
		std::string fileName = oNode->GetAttribute( "file" );
		std::cout<<"  - Loading the Gamma Width image (DIAMOND) from <" << fileName.c_str() << ">"<<std::endl;
		typedef itk::ImageFileReader<KappasImageType> GammaWidthImageReaderType;
		GammaWidthImageReaderType::Pointer reader = GammaWidthImageReaderType::New();
		reader->SetFileName( fileName );
		reader->Update();
		
		this->m_GammaWidthImage = reader->GetOutput();
	}
 

	void MFMIO::WriteHeader( const std::string& fileName )
	{
		crl::FileName baseFileName(fileName);

		//-----------------------------------------
		// Creates the root node
		//-----------------------------------------
		itk::DOMNode::Pointer docDOM = itk::DOMNode::New();
		docDOM->SetName( "MFM" );

		//-----------------------------------------
		// Adds the B=0 image
		//-----------------------------------------
		if ( m_B0Image.IsNotNull() )
			docDOM->AddChildAtEnd( createXmlNode_SingleFile("B0", baseFileName.getCompleteFilePath_WithSufix( "_b0" ))  );

		//-----------------------------------------
		// Creates the tensor list
		//-----------------------------------------
		std::vector<std::string> tensorFiles;
		for ( unsigned int i=0; i<m_TensorImages.size(); i++ )
		{
			itk::FancyString suffix; suffix << "_t" << i;
			std::string ofilename = baseFileName.getCompleteFilePath_WithSufix(suffix);
			tensorFiles.push_back(ofilename);
		}
		docDOM->AddChildAtEnd( createXmlNode_Tensors(tensorFiles) );

		//-----------------------------------------
		// Adds the fractions
		//-----------------------------------------
		if ( m_FractionsImage.IsNotNull() )
			docDOM->AddChildAtEnd( createXmlNode_SingleFile("Fractions", baseFileName.getCompleteFilePath_WithSufix( "_fractions" ))  );

		//-----------------------------------------
		// Adds the DIAMOND Gamma Width image
		//-----------------------------------------
		if ( m_GammaWidthImage.IsNotNull() )
			docDOM->AddChildAtEnd( createXmlNode_SingleFile("DIAMONDGamma", baseFileName.getCompleteFilePath_WithSufix( "_gamma" ))  );

		//-----------------------------------------
		// Adds the D_iso image
		//-----------------------------------------
		if ( m_DisoImage.IsNotNull() )
			docDOM->AddChildAtEnd( createXmlNode_SingleFile("Diso", baseFileName.getCompleteFilePath_WithSufix( "_unrestricted" ))  );

		//-----------------------------------------
		// Now writes the XML file!
		//-----------------------------------------
		baseFileName.setExtension("xml");
		itk::DOMNodeXMLWriter::Pointer xmlWriter = itk::DOMNodeXMLWriter::New();
		xmlWriter->SetInput( docDOM );
		xmlWriter->SetFileName( baseFileName.getCompleteFilePath() );
		xmlWriter->Update();	
	}


	itk::DOMNode::Pointer MFMIO::createXmlNode_Tensors(const std::vector<std::string> & filenames ) 
	{
		itk::DOMNode::Pointer node = itk::DOMNode::New();
		node->SetName( "Tensors" );

		for ( unsigned int i=0; i<filenames.size(); i++ )
		{
			//-----------------------------------------
			// One new node per tensor file
			//-----------------------------------------
			itk::DOMNode::Pointer tensorFile = itk::DOMNode::New();
			tensorFile->SetName( "Tensor" );
			tensorFile->SetAttribute( "file", filenames[i] );
			node->AddChildAtEnd( tensorFile );
		}
		return node;
	}

	itk::DOMNode::Pointer MFMIO::createXmlNode_SingleFile( const std::string& tagname, const std::string& filename ) 
	{
		itk::DOMNode::Pointer node = itk::DOMNode::New();
		node->SetName( tagname );
		node->SetAttribute( "file", filename );
		return node;
	}



	// ESO

    void MFMIO::reshapeLike(MFMIO& MFMIOtemplate)
    {
        const unsigned int nbTensor = MFMIOtemplate.getNumberOfTensors();


        this->m_TensorImages.resize(nbTensor);
        for(unsigned int i = 0; i < nbTensor; i++)
        {
            this->m_TensorImages[i] = TensorImageType::New();
            this->m_TensorImages[i]->SetRegions(MFMIOtemplate.getTensorImage(i)->GetLargestPossibleRegion());
            this->m_TensorImages[i]->SetOrigin(MFMIOtemplate.getTensorImage(i)->GetOrigin());
            this->m_TensorImages[i]->SetSpacing(MFMIOtemplate.getTensorImage(i)->GetSpacing());
            this->m_TensorImages[i]->SetDirection(MFMIOtemplate.getTensorImage(i)->GetDirection());
            this->m_TensorImages[i]->Allocate();

            TensorImageType::PixelType nullT;
            nullT.Fill(0);
            for ( TensorImageRegionIterator it(this->m_TensorImages[i], this->m_TensorImages[i]->GetLargestPossibleRegion()) ; !it.IsAtEnd(); ++it ){
                it.Set(nullT);
            }
        }

        if (MFMIOtemplate.getFractionsImage().IsNotNull())
        {
            this->m_FractionsImage = FractionsImageType::New();
            this->m_FractionsImage->SetRegions(MFMIOtemplate.getFractionsImage()->GetLargestPossibleRegion());
            this->m_FractionsImage->SetOrigin(MFMIOtemplate.getFractionsImage()->GetOrigin());
            this->m_FractionsImage->SetSpacing(MFMIOtemplate.getFractionsImage()->GetSpacing());
            this->m_FractionsImage->SetDirection(MFMIOtemplate.getFractionsImage()->GetDirection());
            this->m_FractionsImage->SetNumberOfComponentsPerPixel(nbTensor);
            this->m_FractionsImage->Allocate();

            FractionsImageType::PixelType initValue(nbTensor);
            initValue.Fill( 0.0 );
            for ( FractionsImageRegionIterator it(this->m_FractionsImage, this->m_FractionsImage->GetLargestPossibleRegion()) ; !it.IsAtEnd(); ++it ){
                it.Set(initValue);
            }
        }

        //todo other ?

    }


    // Save functions (file writer, and name suffix)
    void MFMIO::saveAll( const std::string& filePrefix)
    {

        if(this->getNumberOfTensors() > 0)
        {
            this->saveTensors(filePrefix);
        }

        if ( this->m_FractionsImage.IsNotNull())
        {
            std::string fractionFileName = filePrefix + this->fractionSuffix();
            this->saveFractions(fractionFileName);
        }

        //todo other ?

    }

    void MFMIO::saveTensors(const std::string& filePrefix)
    {
        const unsigned int nbTensor = this->getNumberOfTensors();

        this->m_TensorImages.resize(nbTensor);
        for(unsigned int i = 0; i < nbTensor; i++)
        {
            typedef typename itk::ImageFileWriter<TensorImageType > TensorWriterType;
            TensorWriterType::Pointer tensorWriter = TensorWriterType::New();

            std::string tensorFileName = filePrefix + this->tensorSuffix(i);

            tensorWriter->SetFileName( tensorFileName );
            tensorWriter->SetInput( this->m_TensorImages[i] );
            tensorWriter->SetUseCompression(true); // todo parameter ?
            tensorWriter->Update();
        }

    }

    void MFMIO::saveFractions( const std::string& fractionName )
    {

        typedef typename itk::ImageFileWriter<FractionsImageType > FractionsWriterType;
        FractionsWriterType::Pointer fractionWriter = FractionsWriterType::New();

        fractionWriter->SetFileName( fractionName );
        fractionWriter->SetInput( this->m_FractionsImage );
        fractionWriter->SetUseCompression(true); // todo parameter
        fractionWriter->Update();
    }


    std::string MFMIO::fractionSuffix()
    {
        return std::string("_fraction.nrrd");
    }

    std::string MFMIO::tensorSuffix(unsigned int i)
    {
        return (std::string("_t") + std::to_string(i) + std::string(".nrrd"));
    }

}
