/**********************************************************************************************//**
 * \file	crlMFMIO.h
 *
 * \brief	Declares the crl::MFMIO class. 
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


#ifndef h_CRL_MFMIO
#define h_CRL_MFMIO

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkDiffusionTensor3D.h>
#include "itkImageRegion.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <itkDOMNode.h>

namespace crl {

	/**********************************************************************************************//**
	 * \class	MFMIO
	 *
	 * \brief	
	 *
	 * \author	Benoit Scherrer
	 * \date	August 2013
	 **************************************************************************************************/
	class MFMIO
	{
	public:
		// Types definition
		typedef itk::DiffusionTensor3D< float >				TensorPixelType;
		typedef itk::Image< TensorPixelType, 3 >			TensorImageType;
		typedef itk::VectorImage<float,3>					FractionsImageType;
		typedef itk::VectorImage<float,3>					KappasImageType;
		typedef itk::Image<float,3>							FloatImageType;

		typedef std::vector<TensorImageType::Pointer>		TensorImageListType;


	    typedef itk::ImageRegionIterator< TensorImageType >         TensorImageRegionIterator;
	    typedef itk::ImageRegionIterator< FractionsImageType >      FractionsImageRegionIterator;

	public:
		MFMIO();
		~MFMIO();

		int					getNumberOfTensors() const;
		bool				freeDiffusion() const;

		const TensorImageListType&						getTensorImages();
		const TensorImageType::Pointer&					getTensorImage(int id) ;

		const FractionsImageType::Pointer&				getFractionsImage() ;
		const KappasImageType::Pointer&					getKappasImage() ;
		const FloatImageType::Pointer&					getB0Image() ;
		const FloatImageType::Pointer&					getDisoImage() ;


		void	ReadFromCommandLine( const std::vector<std::string>& fileNames, const std::string& fractionName="", const std::string& kappaName="" );


		void	ReadXML( const std::string& fileName );
		void	WriteHeader( const std::string& fileName );

		static itk::DOMNode::Pointer createXmlNode_Tensors(const std::vector<std::string> & filenames ) ;
		static itk::DOMNode::Pointer createXmlNode_SingleFile( const std::string& tagname, const std::string& filename ) ;

		void	clearTensorImages();
		void	addTensorImage(const TensorImageType::Pointer& image );

		void	setFractionsImage(const FractionsImageType::Pointer& image );
		void	setKappasImage(const KappasImageType::Pointer& image );
		void	setB0Image(const FloatImageType::Pointer& image );
		void	setDisoImage(const FloatImageType::Pointer& image );


	protected:
		void	loadTensors( itk::DOMNode::Pointer node );
		void	loadFractions( itk::DOMNode::Pointer node );
		void	loadB0( itk::DOMNode::Pointer node );
		void	loadDiso( itk::DOMNode::Pointer node );
		void	loadDIAMOND( itk::DOMNode::Pointer node );

	protected:
		TensorImageListType							m_TensorImages ;
		FractionsImageType::Pointer					m_FractionsImage;
		KappasImageType::Pointer				    m_GammaWidthImage;
		FloatImageType::Pointer						m_B0Image;
		FloatImageType::Pointer						m_DisoImage;

	public:
		void reshapeLike(MFMIO& MFMIOtemplate);

        void saveAll( const std::string& filePrefix);
        void saveTensors( const std::string& filePrefix);
        void saveFractions( const std::string& fractionName );

	private:
        std::string fractionSuffix();
        std::string tensorSuffix(unsigned int i);
	};

	

}

#endif
