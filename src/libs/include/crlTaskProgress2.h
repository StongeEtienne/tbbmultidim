/**********************************************************************************************//**
 * \file	crlTaskProgress.h
 *
 * \brief	Declares the crl::TaskProgress class. 
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

#ifndef h_CRL_TASK_PROGRESS2
#define h_CRL_TASK_PROGRESS2

#include <string>
#include <itkTimeProbe.h>
#include <itkCommand.h>
#include <itkProcessObject.h>

namespace crl {

	
	class TaskProgress2
	{
	public:
		TaskProgress2( unsigned long int counterMaxVal=0, bool showRemainingTime=true );
		~TaskProgress2();
		
		void		InitializeProgress(unsigned long int counterMaxVal=0, bool startTimer=true);
		void		StartProgress();
		void		IncrementProgress();
		
		void		SetShowEstimatedRemainingTime(bool b);
		bool		GetShowEstimatedRemainingTime() const;

		static		std::string ConvertTimeToString(itk::RealTimeClock::TimeStampType t) ;

		void		SetLinePrefix(const std::string& prefix);
		std::string GetLinePrefix() const;
		void		SetLineSuffix(const std::string& suffix);
		std::string GetLineSuffix() const;

		void		SetSilent(bool b) ;
		bool		GetSilent() const;

		void		StdOutLock();
		void		StdOutUnlock();

	protected:
		void		ShowProgress();

	private:
		class		PrivateData;
		PrivateData *d;

	};

	/**********************************************************************************************//**
	 * \class	ItkFilterProgressObserver
	 *
	 * \brief	Defines a itk::Command that can be plugged to a ITK filter to show the progress of
	 * 			the filter.
	 * 			
	 * 			EXAMPLE:  
	 * 			ItkFilterProgressObserver::Pointer progressObs = ItkFilterProgressObserver::New();  
	 * 			progressObs->SetProgressLinePrefix("Processing... ");  
	 *			filter->AddObserver(ProgressEvent(), progressObs );  
	 * 			filter->AddObserver(StartEvent(), progressObs );  
	 *	
	 * \author	Benoit Scherrer
	 * \date	June 2010
	*************************************************************************************************/
	//class ItkFilterProgressObserver2 : public itk::Command 
	//{
	//public:
	//	typedef  ItkFilterProgressObserver2  Self;
	//	typedef  itk::Command				Superclass;
	//	typedef  itk::SmartPointer<Self>	Pointer;
	//	typedef   itk::ProcessObject		ProcessType;
	//	typedef   const ProcessType   *		ProcessPointer;

	//	itkNewMacro( Self );

	//public:
	//	void	Execute(itk::Object *caller, const itk::EventObject & event);
	//	void	Execute(const itk::Object * object, const itk::EventObject & event);
	//	void	SetProgressLinePrefix(const std::string& str);

	//protected:
	//	ItkFilterProgressObserver();
	//	TaskProgress m_TaskProgress;
	//};

}

#endif
