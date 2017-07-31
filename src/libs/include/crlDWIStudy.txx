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


#ifndef CRL_DWI_STUDY_TXX
#define CRL_DWI_STUDY_TXX

#include "crlFileName.cxx"
#include "crlFileUtil.cxx"

#include "crlDWIStudy.h"
#include "crlNHDRFileUtils.h"
#include "crlFileUtil.h"
#include "crlRobustMeanImageFilter.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientImageFilter.h>

#include <iostream>
#ifdef WIN32
#include <cctype>
#endif

namespace crl{

/**********************************************************************************************//**
 * \fn	template <class TPixel> DWIStudy<TPixel>::DWIStudy()
 *
 * \brief	Constructor
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type. 
*************************************************************************************************/
template <class TPixel>
DWIStudy<TPixel>::DWIStudy():
m_GradientCoordinatesMode(GRADIENTS_SCANNER_COORDINATES)
{
	NominalBValue = -1;
}

/**********************************************************************************************//**
 * \fn	DWIStudy::~DWIStudy()
 *
 * \brief	Destructor. 
 *
 * \author	Benoit Scherrer
 * \date	February 2011
*************************************************************************************************/
template <class TPixel>
DWIStudy<TPixel>::~DWIStudy()
{
	Gradients.clear();
	BValues.clear();
	OriginalGradients.clear();
}
		
template <class TPixel>
DWIStudy<TPixel>::DWIStudy(const DWIStudy<TPixel>&ref)
{
	this->fileName						= ref.fileName;
	this->DWIData						= ref.DWIData;
	this->Gradients						= ref.Gradients;
	this->BValues						= ref.BValues;
	this->NominalBValue					= ref.NominalBValue;
	this->OriginalGradients				= ref.OriginalGradients;	
	this->m_GradientCoordinatesMode		= ref.m_GradientCoordinatesMode;
}
		
template <class TPixel>
void DWIStudy<TPixel>::operator=(const DWIStudy<TPixel>& ref)
{
	this->fileName						= ref.fileName;
	this->DWIData						= ref.DWIData;
	this->Gradients						= ref.Gradients;
	this->BValues						= ref.BValues;
	this->NominalBValue					= ref.NominalBValue;
	this->OriginalGradients				= ref.OriginalGradients;	
	this->m_GradientCoordinatesMode		= ref.m_GradientCoordinatesMode;

}

template <class TPixel>	
void DWIStudy<TPixel>::Allocate( const DWIStudy<TPixel>& geometry )
{
	this->Allocate(geometry.DWIData );
}
	
template <class TPixel>
void DWIStudy<TPixel>::Allocate( const DWIStudy<TPixel>& geometry, int numberOfGradients )
{
	this->Allocate(geometry.DWIData, numberOfGradients );
}

template <class TPixel>
void DWIStudy<TPixel>::Allocate( const typename DWIGradientImageSetType::Pointer& geometry )
{
	this->DWIData = DWIGradientImageSetType::New();
	this->DWIData->SetRegions(geometry->GetLargestPossibleRegion());
	this->DWIData->SetOrigin(geometry->GetOrigin());
	this->DWIData->SetSpacing(geometry->GetSpacing());
	this->DWIData->SetDirection(geometry->GetDirection());
	this->DWIData->SetNumberOfComponentsPerPixel(geometry->GetNumberOfComponentsPerPixel()); 
	this->DWIData->Allocate();
	this->fileName="";
}

template <class TPixel>
void DWIStudy<TPixel>::Allocate( const typename DWIGradientImageSetType::Pointer& geometry, int numberOfGradients )
{
	this->DWIData = DWIGradientImageSetType::New();
	this->DWIData->SetRegions(geometry->GetLargestPossibleRegion());
	this->DWIData->SetOrigin(geometry->GetOrigin());
	this->DWIData->SetSpacing(geometry->GetSpacing());
	this->DWIData->SetDirection(geometry->GetDirection());
	this->DWIData->SetNumberOfComponentsPerPixel(numberOfGradients); 
	this->DWIData->Allocate();
	this->fileName="";
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::LoadStudy ( const std::string& fileName )
 *
 * \brief	Loads a DWI study. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type.  
 * \param	fileName	Filename of the file. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::LoadStudy ( const std::string& fileName )
{
	crl::FileName file(fileName);
	std::string fileExt=file.getExtension();

	if ( fileExt=="nhdr" )
		LoadNHDRStudy(fileName);
	else if ( fileExt=="nii" || fileExt=="nii.gz" || fileExt=="hdr" )
		LoadFSLStudy(fileName);
	else {
		std::string errMsg = "Error while opening "+fileName+".\nDon't know how to open a '"+fileExt+"' DWI study"; 
		throw itk::ExceptionObject(__FILE__,__LINE__,errMsg, "");
	}

	this->fileName=fileName;
}

template <class TPixel>
void DWIStudy<TPixel>::LoadFSLStudy ( const std::string& fileName )
{
	//---------------------------------------
	// First load the DWI vector image
	//---------------------------------------
	typedef itk::Image< TPixel, 4 > Image4DType;
	std::cout<<"- Loading the DWI file <" << fileName << ">..." << std::endl;
	typedef itk::ImageFileReader<Image4DType> Image4DReaderType;
	typename Image4DReaderType::Pointer reader = Image4DReaderType::New();
	reader->SetFileName(fileName);
	reader->Update();

	this->DWIData = crl::ConvertNDToVectorImage<TPixel,4>(reader->GetOutput());

	//---------------------------------------
	// Now read vectors from either
	//    [filename].bvecs
	//    [filename].bvec
	//    bvecs
	//    bvec
	//---------------------------------------
	crl::FileName filename(fileName);
	std::string bvecFileName;
	if ( crl::file_exists(filename.getWithNewExtension("bvecs")))
		bvecFileName = filename.getWithNewExtension("bvecs");
	else if ( crl::file_exists(filename.getWithNewExtension("bvec")))
		bvecFileName = filename.getWithNewExtension("bvec");
	else if ( crl::file_exists(filename.getWithNewFilenameAndExt("bvec")))
		bvecFileName = filename.getWithNewFilenameAndExt("bvec");
	else
		bvecFileName = filename.getWithNewFilenameAndExt("bvecs");

	std::cout<<"- Read gradient vectors from <"<<bvecFileName.c_str()<<">..."<<std::endl;
	std::ifstream ifile_bvecs;
	ifile_bvecs.open ( bvecFileName.c_str(), ifstream::in);
	if  ( !ifile_bvecs.is_open() ) 
		throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot read/find the gradient vector file for reading", "");
	
	std::string sgX, sgY, sgZ;
	std::getline(ifile_bvecs, sgX);
	std::getline(ifile_bvecs, sgY);
	std::getline(ifile_bvecs, sgZ);

	std::vector<double> vgX, vgY, vgZ;
	string2vector(sgX, vgX );
	string2vector(sgY, vgY );
	string2vector(sgZ, vgZ );
	ifile_bvecs.close();


	// Check the different sizes
	if ( vgX.size()!=vgY.size() ||  vgX.size()!=vgZ.size() || vgY.size()!=vgZ.size())
	{
		std::cout<<"gX (" << vgX.size() <<") :" ;
		for ( unsigned int i=0; i<vgX.size(); i++ )
			std::cout<<" " << vgX[i];

		std::cout<< std::endl<<"gY (" << vgY.size() <<") :" ;
		for ( unsigned int i=0; i<vgY.size(); i++ )
			std::cout<<" " << vgY[i];

		std::cout<< std::endl<<"gZ (" << vgZ.size() <<") :" ;
		for ( unsigned int i=0; i<vgZ.size(); i++ )
			std::cout<<" " << vgZ[i];
		std::cout<< std::endl;

		throw itk::ExceptionObject(__FILE__,__LINE__,"Error while reading the gradient vectors. The X,Y and Z components doesn't have the same size.", "");
	}		

	// If invalid, try to open with transposition
	if ( vgX.size()==0 || vgY.size()==0 || vgZ.size()==0 )
	{
                std::cout<<"  Error. Try to read the vectors line by line." << std::endl;
		vgX.clear();
		vgY.clear(); 
		vgZ.clear();
		ifile_bvecs.open ( bvecFileName.c_str(), ifstream::in);
		if  ( !ifile_bvecs.is_open() ) 
			throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot read/find the gradient vector file for reading", "");

		std::string line;
		while (std::getline(ifile_bvecs, line))
		{
			std::istringstream iss(line);
			double x, y, z;

			if (!(iss >> x >> y >> z )) { break; } // error

			vgX.push_back(x);
			vgY.push_back(y);
			vgZ.push_back(z);
			// process pair (a,b)
		}

		ifile_bvecs.close();
	}

	std::cout<<"  (found "<<vgX.size()<<" gradient vectors)"<<std::endl;

	// Construct a std::vector<GradientDirectionType>
	this->Gradients.clear();
	for ( unsigned int i=0; i<vgX.size() ; i++ )
	{	
		GradientVectorType v;
		v[0] = vgX[i];
		v[1] = vgY[i];
		v[2] = vgZ[i];
		this->Gradients.push_back(v);
	}

	//-------------------------------------
	//  Read the B-Values from either
	//    [filename].bvals
	//    [filename].bvsl
	//    bvals
	//    bval
	//-------------------------------------
	std::string bvalFileName;
	if ( crl::file_exists(filename.getWithNewExtension("bvals")))
		bvalFileName = filename.getWithNewExtension("bvals");
	else if ( crl::file_exists(filename.getWithNewExtension("bval")))
		bvalFileName = filename.getWithNewExtension("bval");
	else if ( crl::file_exists(filename.getWithNewFilenameAndExt("bval")))
		bvalFileName = filename.getWithNewFilenameAndExt("bval");
	else
		bvalFileName = filename.getWithNewFilenameAndExt("bvals");

	cout<<"- Read B-Values from <"<<bvalFileName.c_str()<<">..."<<endl;
	std::ifstream ifile_bvals;
	ifile_bvals.open ( bvalFileName.c_str(), ifstream::in);
	if  ( !ifile_bvals.is_open() ) 
	{
		throw itk::ExceptionObject(__FILE__,__LINE__,"Cannot read/find the b-values file for reading", "");
	}

	std::string sgB;
	std::getline(ifile_bvals, sgB);

	this->BValues.clear();
	string2vector(sgB, this->BValues );

	ifile_bvals.close();

	std::cout<<"  (found "<<this->BValues.size()<<" b-values)"<<std::endl;
	 
	// Check the size
	if ( this->BValues.size()!=vgX.size()  )
	{
		std::cout<<"B-Values: (" << this->BValues.size() <<") :"  ;
		for ( unsigned int i=0; i<this->BValues.size(); i++ )
			std::cout<<" " << this->BValues[i] ;

		std::cout<<std::endl<<"Gradients: (" << vgX.size() <<") :" <<std::endl ;
		for ( unsigned int i=0; i<vgX.size(); i++ )
			std::cout<<"  [" << vgX[i]<<","<<vgY[i]<<","<<vgZ[i]<<"]"<<std::endl ;;
		std::cout<<std::endl;

		throw itk::ExceptionObject(__FILE__,__LINE__,"Error while reading the B-values. The number of b-values is not equal to the number of gradient directions.", "");
	}
}


template <class TPixel>
bool DWIStudy<TPixel>::is_number(const std::string& str)
{
	std::string::const_iterator it = str.begin();
	while ( it!=str.end() && (std::isdigit(*it) || *it=='.' || *it=='-' || *it=='+' || *it=='e') ) ++it;
	return !str.empty() && it==str.end(); 
}

	/**********************************************************************************************//**
	 * \fn	void string2vector(const std::string& str, vector<double> &outputVals )
	 *
	 * \brief	Converts a space-separated list of numbers in a string to a vector<double> 
	 *
	 * \param	str					The string. 
	 * \param [in,out]	outputVals	The output vals. 
	*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::string2vector(const std::string& str, std::vector<double> &outputVals )
{
	stringstream ss(str+std::string(" \n"));

	string val;
	outputVals.clear();

	getline(ss,val,' ');
	while ( ss.rdstate() == ios::goodbit )
	{
		if ( is_number(val) )
			outputVals.push_back( atof(val.c_str()) );
		getline(ss,val,' ');
	}
	if ( val!="" && val!=" " )
	{
		if ( is_number(val) )
			outputVals.push_back( atof(val.c_str()) );		//don't forget the last one
	}

	int value;
	ss>>value;
}

template <class TPixel>
void DWIStudy<TPixel>::LoadNHDRStudy ( const std::string& fileName )
{
	//---------------------------------------
	// First load the DWI vector image
	//---------------------------------------
	std::cout<<"- Loading the DWI file <" << fileName << ">..." << std::endl;
	typedef itk::ImageFileReader<DWIGradientImageSetType> DWIStudyReaderType;
	typename DWIStudyReaderType::Pointer reader = DWIStudyReaderType::New();
	reader->SetFileName(fileName);
	reader->Update();
	this->DWIData = reader->GetOutput();

	//--------------------------------------------
	// Prepare to read the B-Values/B-Vectors
	//--------------------------------------------
	NominalBValue = 0;
	Gradients.clear();
	BValues.clear();

	GradientVectorType vect3d;
	bool m_ReadB0 = false;
	int nbNullGrad = 0;

	//--------------------------------------------
	// Read the keys in the meta dictionary
	//--------------------------------------------
	std::vector<std::string> imgMetaKeys = DWIData->GetMetaDataDictionary().GetKeys();
	std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
	std::string metaString;
	for (; itKey != imgMetaKeys.end(); itKey ++)
	{
		itk::ExposeMetaData<std::string> (DWIData->GetMetaDataDictionary(), *itKey, metaString);

		//--------------------------------------------
		// Is it a gradient vector?
		//--------------------------------------------
		if (itKey->find("DWMRI_gradient") != std::string::npos)
		{ 
			double x,y,z;
			sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
			vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;

			if ( vect3d.one_norm()<1e-8 ) nbNullGrad++;
			Gradients.push_back(vect3d);
		}

		//--------------------------------------------
		// Is it the nominal B-value?
		//--------------------------------------------
		else if (itKey->find("DWMRI_b-value") != std::string::npos)
		{
			m_ReadB0 = true;
			NominalBValue = atof(metaString.c_str());
		}
	}

	//--------------------------------------------
	// Check that the nominal b-value has been read
	//--------------------------------------------
	if(!NominalBValue)
		throw itk::ExceptionObject(__FILE__,__LINE__,"Nominal B-Value not found in the NHDR file", "");

	std::cout << "  - Found " << Gradients.size()-nbNullGrad << " gradient vectors, "<<nbNullGrad<<" B=0 and a nominal B-value of "<<NominalBValue<<"."<<std::endl;

	//--------------------------------------------
	// Now normalize the gradients and compute the b-values
	//--------------------------------------------
	for ( unsigned int i=0; i<Gradients.size(); i++ )
	{
		vect3d = Gradients[i];

		if ( vect3d.one_norm() <= 1e-8 )
		{
			vect3d[0]=vect3d[1]=vect3d[2]=0;
			Gradients[i] = vect3d;
			BValues.push_back(0);
		}
		else
		{
			double nn = vect3d.two_norm();
			double b = NominalBValue * nn * nn;
			Gradients[i] = ( vect3d / nn );
			BValues.push_back(b);
		}
	}

	//----------------------------------------
	// Keep the original gradients in case OrientStudyInAxial is called
	//----------------------------------------
	OriginalGradients.clear();
	for ( unsigned int i=0; i<Gradients.size(); i++ )
		OriginalGradients.push_back(Gradients[i]);
}

template <class TPixel>
void DWIStudy<TPixel>::Set(const typename DWIGradientImageSetType::Pointer& data, double nominalB, const GradientVectorSetType& normalizedGradients, const BValueSetType& bValues )
{
	this->DWIData = data;
	this->NominalBValue = nominalB;
	this->Gradients = normalizedGradients;
	this->BValues = bValues;
	this->OriginalGradients = normalizedGradients;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::WriteStudy( const std::string& fileName )
 *
 * \brief	Writes a DWI study. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	Pixel type.  
 * \param	fileName	Filename of the file. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::WriteStudy( const std::string& fileName, bool useCompression )
{

	crl::NHDR_WriteFromVectorImage<TPixel, 3> ( DWIData,
		Gradients,
		BValues,
		fileName,
		useCompression
		);

	this->fileName = fileName;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::WriteStudyWithNewData( const std::string& fileName,
 * 		typename DWIGradientImageSetType::Pointer data )
 *
 * \brief	Writes the study with new data but by keeping the same gradient vectors/b-values. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 * \param	fileName	Filename of the file. 
 * \param	data		The data. 
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::WriteStudyWithNewData( const std::string& fileName, typename DWIGradientImageSetType::Pointer data, bool useCompression )
{
	crl::NHDR_WriteFromVectorImage<TPixel, 3> ( data,
		Gradients,
		BValues,
		fileName,
		useCompression
		);

	this->fileName = fileName;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::SetGradientCoordinates(GradientCoordinatesType gc)
 *
 * \brief	Sets the gradient coordinates mode. If GRADIENTS_MFRAME_COORDINATES, the
 * 			OrientStudyInAxial function will reorient the gradients to the scanner coordinates. 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \typeparam	TPixel	. 
 * \param	gc	The gradient coordinates type. 
 *
 * \sa	OrientStudyInAxial
*************************************************************************************************/
template <class TPixel>
void 
DWIStudy<TPixel>::SetGradientCoordinates(GradientCoordinatesType gc)
{
	m_GradientCoordinatesMode = gc;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> typename DWIStudy<TPixel>::GradientCoordinatesType DWIStudy<TPixel>::GetGradientCoordinates()
 *
 * \brief	Gets the current gradient coordinates. mode 
 *
 * \author	Benoit Scherrer
 * \date	December 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The gradient coordinates mode. 
*************************************************************************************************/
template <class TPixel>
typename DWIStudy<TPixel>::GradientCoordinatesType
DWIStudy<TPixel>::GetGradientCoordinates()
{
	return m_GradientCoordinatesMode;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetNumberOfB0() const
 *
 * \brief	Gets the number of b=0 images
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of b=0 images. 
*************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetNumberOfB0() const
{
	int nb=0;
	for ( unsigned int i=0; i<BValues.size(); i++ )
	{
		if ( BValues[i]==0 ) nb++;
	}
	return nb;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetNumberOfNonNullVectors() const
 *
 * \brief	Gets the number of non null gradient vectors. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \return	The number of non null gradient vectors. 
*************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetNumberOfNonNullVectors() const
{
	int nb=0;
	for ( unsigned int i=0; i<BValues.size(); i++ )
	{
		if ( BValues[i]!=0 ) nb++;
	}
	return nb;
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> unsigned int DWIStudy<TPixel>::GetTotalNumberOfImages() const
 *
 * \brief	Gets the total number of images.
 *
 * \author	Benoit Scherrer
 * \date	June 2011
 *
 * \tparam	TPixel	Type of the pixel.
 *
 * \return	The total number of images.
 **************************************************************************************************/
template <class TPixel>
unsigned int DWIStudy<TPixel>::GetTotalNumberOfImages() const
{
	return BValues.size();
}

/**********************************************************************************************//**
 * \fn	template <class TPixel> double DWIStudy<TPixel>::GetMaximumBValue() const
 *
 * \brief	Gets the maximum b value.
 *
 * \author	Benoit Scherrer
 * \date	June 2011
 *
 * \exception	itk::ExceptionObject	Thrown if error
 *
 * \tparam	TPixel	Type of the pixel.
 *
 * \return	The maximum b value.
 **************************************************************************************************/
template <class TPixel>
double DWIStudy<TPixel>::GetMaximumBValue() const
{
	if (BValues.size()<1 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. The array of b-values is empty.", "DWIStudy::GetMaximumBValue");

	double max = BValues[0];
	for ( unsigned int i=0; i<BValues.size(); i++ )
		if (BValues[i]>max ) max = BValues[i];

	return max;
}

template <class TPixel>
void DWIStudy<TPixel>::ComputeScaledGradients( GradientVectorSetType& scaledGradients ) const
{
	if (BValues.size()!=Gradients.size())
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. The number of b-values does not match the number of gradients", "");

	double maxB = GetMaximumBValue();

	scaledGradients.clear();
	GradientVectorType v;
	for ( unsigned int i=0; i<Gradients.size(); i++ )
	{
		v = Gradients[i] * sqrt(BValues[i]/maxB);
		scaledGradients.push_back(v);
	}
}


/**********************************************************************************************//**
 * \fn	template <class TPixel> void DWIStudy<TPixel>::OrientStudyInAxial ()
 *
 * \brief	Orients the current study in axial. If the gradient coordinates mode is set to
 * 			GRADIENTS_MFRAME_COORDINATES, the gradient vectors are modified accordingly to the
 * 			axes flip and permutations performed by the orienter. 
 *
 * \author	Benoit Scherrer
 * \date	September 2010
 *
 * \typeparam	TPixel	. 
 *
 * \sa	SetGradientCoordinates
*************************************************************************************************/
template <class TPixel>
void DWIStudy<TPixel>::OrientStudyInAxial ()
{
	//-----------------------------------------
	// Orient the input image
	//-----------------------------------------
	std::cout<<"  - Orient the study to axial..."<<std::endl;
	typedef typename itk::OrientImageFilter<DWIGradientImageSetType, DWIGradientImageSetType> OrientVectorImageFilterType;
	typename OrientVectorImageFilterType::Pointer orienter = OrientVectorImageFilterType::New();
	orienter->SetInput(DWIData);
	orienter->UseImageDirectionOn();
	orienter->SetDesiredCoordinateOrientationToAxial();
	orienter->Update();
	DWIData = orienter->GetOutput();
	orienter->GetOutput()->DisconnectPipeline();
	DWIData->Update();

	//-------------------------------------------------------
	// Apply the permutation and then the flip performed by the orienter to the vectors
	//-------------------------------------------------------
	if ( m_GradientCoordinatesMode==GRADIENTS_MFRAME_COORDINATES )
	{
		std::cout<<"  - Apply the permutations/flips to the gradients..."<<std::endl;
		const typename OrientVectorImageFilterType::PermuteOrderArrayType & permuteOrder = orienter->GetPermuteOrder();
		const typename OrientVectorImageFilterType::FlipAxesArrayType & flippedAxes = orienter->GetFlipAxes();
		GradientVectorType outVector;
		for ( unsigned int v=0; v<Gradients.size(); v++ )
		{
			// Permutations
			for ( unsigned int i=0;i<3 ; i++ )
			{
				GradientVectorType inVector = Gradients[v];
				outVector[i] = inVector[permuteOrder[i]];
			}
			// Flips
			for ( unsigned int i=0;i<3 ; i++ )
			{
				if ( flippedAxes[i] && outVector[i]!=0 ) outVector[i] = -outVector[i];
			}
			Gradients[v] = outVector;
		}
	}
}

template <class TPixel>
typename DWIStudy<TPixel>::DWIGradientImageType::Pointer DWIStudy<TPixel>::ComputeMeanB0()
{
	//-----------------------------------------
	// Extract images
	//-----------------------------------------
	std::cout<<"- Extract the B=0 images..."<<std::endl;
	typedef crl::RobustMeanImageFilter< DWIGradientImageType, DWIGradientImageType >  RobustMeanImageFilterType;
	typename RobustMeanImageFilterType::Pointer meanFilter = RobustMeanImageFilterType::New();

	//-------------------------------------------------------
	// Add all B!=0 images and connect them to the robust mean filter
	//-------------------------------------------------------
	int numImg = 0; 
	for( unsigned int b=0; b<BValues.size(); b++ )
	{
		// Take only zero-norm gradients
		if( BValues[b]==0 ) 
		{
			typename DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3 >( DWIData, b );
			meanFilter->SetInput(numImg, oneImage);
			numImg++;
		}
	}

	//-------------------------------------------------------
	// Checking...
	//-------------------------------------------------------
	if ( numImg==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. No null b-value found", "");

	//-------------------------------------------------------
	// Computes!
	//-------------------------------------------------------
	meanFilter->Update();
	std::cout<<std::endl;
	return meanFilter->GetOutput();
}


template <class TPixel>
typename DWIStudy<TPixel>::DWIGradientImageType::Pointer DWIStudy<TPixel>::ComputeMeanBNonNull()
{
	//-----------------------------------------
	// Extract images
	//-----------------------------------------
	std::cout<<"- Extract the B!=0 images..."<<std::endl;
	typedef crl::RobustMeanImageFilter< DWIGradientImageType, DWIGradientImageType >  RobustMeanImageFilterType;
	typename RobustMeanImageFilterType::Pointer meanFilter = RobustMeanImageFilterType::New();
	meanFilter->SetShowTaskProgress(false);

	//-------------------------------------------------------
	// Add all B!=0 images and connect them to the robust mean filter
	//-------------------------------------------------------
	int numImg = 0; 
	for( unsigned int b=0; b<BValues.size(); b++ )
	{
		// Take only zero-norm gradients
		if( BValues[b]!=0 ) 
		{
			typename DWIGradientImageType::Pointer oneImage = crl::ExtractComponentOfVectorImage< TPixel, 3 >( DWIData, b );
			meanFilter->SetInput(numImg, oneImage);
			numImg++;
		}
	}

	//-------------------------------------------------------
	// Checking...
	//-------------------------------------------------------
	if ( numImg==0 )
		throw itk::ExceptionObject(__FILE__,__LINE__,"ERROR. No non-null b-value found", "");

	//-------------------------------------------------------
	// Computes!
	//-------------------------------------------------------
	meanFilter->Update();
	std::cout<<std::endl;
	return meanFilter->GetOutput();
}




/**********************************************************************************************//**
 * \fn	const std::string& DWIStudy<TPixel>::GetFileName() const
 *
 * \brief	Gets the file name.
 *
 * \author	Benoit Scherrer
 * \date	February 2014
 *
 * \tparam	TPixel	Type of the pixel.
 *
 * \return	The file name.
 **************************************************************************************************/
template <class TPixel>
const std::string&	DWIStudy<TPixel>::GetFileName() const
{
	return this->fileName;
}

} /* end namespace crl. */

#endif





  




