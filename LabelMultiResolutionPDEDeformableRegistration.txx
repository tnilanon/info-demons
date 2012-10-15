/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkLabelMultiResolutionPDEDeformableRegistration.txx,v $
Language:  C++
Date:      $Date: 2009-01-26 21:45:51 $
Version:   $Revision: 1.33 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef LabelMultiResolutionPDEDeformableRegistration_txx
#define LabelMultiResolutionPDEDeformableRegistration_txx
#include "LabelMultiResolutionPDEDeformableRegistration.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_math.h"

#include "HistogramField.h"


namespace itk {

	/**
	* Default constructor
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::LabelMultiResolutionPDEDeformableRegistration()
	{

		this->SetNumberOfRequiredInputs(2);

		typename DefaultRegistrationType::Pointer registrator =
			DefaultRegistrationType::New();
		m_RegistrationFilter = static_cast<RegistrationType*>(
			registrator.GetPointer() );

		m_MovingImagePyramid  = MovingImagePyramidType::New();
		m_FixedImagePyramid     = FixedImagePyramidType::New();
		m_LabelImagePyramid   = LabelImagePyramidType::New();
		m_MovingGradImagePyramid = GradImagePyramidType::New();
		//m_FixedGradImagePyramid = GradImagePyramidType::New();

		m_FieldExpander     = FieldExpanderType::New();
		m_InitialDeformationField = NULL;

		m_NumberOfLevels = 3;
		m_NumberOfIterations.resize( m_NumberOfLevels );
		m_FixedImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		m_MovingImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		m_FixedHistogramPyramid=new MultiChannelPyramidType(m_NumberOfLevels);
		m_MovingHistogramPyramid=new MultiChannelPyramidType(m_NumberOfLevels);
		//m_FixedHistogramPyramid->SetNumberOfLevels( m_NumberOfLevels );
		//m_MovingHistogramPyramid->SetNumberOfLevels( m_NumberOfLevels );
		//m_LabelImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		//m_MovingGradImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		// m_FixedGradImagePyramid->SetNumberOfLevels( m_NumberOfLevels );

		unsigned int ilevel;
		for( ilevel = 0; ilevel < m_NumberOfLevels; ilevel++ )
		{
			m_NumberOfIterations[ilevel] = 10;
		}
		m_CurrentLevel = 0;

		m_StopRegistrationFlag = false;

	}


	/*
	* Set the moving image image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::SetMovingImage(
		const MovingImageType * ptr )
	{
		this->ProcessObject::SetNthInput( 2, const_cast< MovingImageType * >( ptr ) );
	}


	/*
	* Get the moving image image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	const typename LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::MovingImageType *
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GetMovingImage(void) const
	{
		return dynamic_cast< const MovingImageType * >
			( this->ProcessObject::GetInput( 2 ) );
	}

	/*
	* Set the label image image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::SetLabelImage(
		const LabelImageType * ptr )
	{
		this->ProcessObject::SetNthInput( 3, const_cast< LabelImageType * >( ptr ) );
	}


	/*
	* Get the label image image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	const typename LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::LabelImageType *
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GetLabelImage(void) const
	{
		return dynamic_cast< const LabelImageType * >
			( this->ProcessObject::GetInput( 3 ) );
	}
	/*
	* Set the MovingGradImage.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::SetMovingGradImage(
		const GradImageType * ptr )
	{
		this->ProcessObject::SetNthInput( 4, const_cast< GradImageType * >( ptr ) );
	}


	/*
	* Get the MovingGradImage.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	const typename LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GradImageType *
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GetMovingGradImage(void) const
	{
		return dynamic_cast< const GradImageType * >
			( this->ProcessObject::GetInput( 4 ) );
	}


	/*
	* Set the fixed image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::SetFixedImage(
		const FixedImageType * ptr )
	{
		this->ProcessObject::SetNthInput( 1, const_cast< FixedImageType * >( ptr ) );
	}


	/*
	* Get the fixed image.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	const typename LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::FixedImageType *
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GetFixedImage(void) const
	{
		return dynamic_cast< const FixedImageType * >
			( this->ProcessObject::GetInput( 1 ) );
	}

	/*
	* 
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	std::vector<SmartPointer<DataObject> >::size_type
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GetNumberOfValidRequiredInputs() const
	{
		typename std::vector<SmartPointer<DataObject> >::size_type num = 0;

		if (this->GetFixedImage())
		{
			num++;
		}

		if (this->GetMovingImage())
		{
			num++;
		}
		

		return num;
	}


	/**
	* Set the number of multi-resolution levels
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::SetNumberOfLevels(
		unsigned int num )
	{
		if( m_NumberOfLevels != num )
		{
			this->Modified();
			m_NumberOfLevels = num;
			m_NumberOfIterations.resize( m_NumberOfLevels );
		}

		if( m_MovingImagePyramid && m_MovingImagePyramid->GetNumberOfLevels() != num )
		{
			m_MovingImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		}
		if( m_FixedImagePyramid && m_FixedImagePyramid->GetNumberOfLevels() != num )
		{
			m_FixedImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
		}  

		if(m_FixedHistogramPyramid)
			m_FixedHistogramPyramid->SetNumberOfLevels( m_NumberOfLevels );
		if(m_MovingHistogramPyramid)
			m_MovingHistogramPyramid->SetNumberOfLevels( m_NumberOfLevels );
	}


	/**
	* Standard PrintSelf method.
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "NumberOfLevels: " << m_NumberOfLevels << std::endl;
		os << indent << "CurrentLevel: " << m_CurrentLevel << std::endl;

		os << indent << "NumberOfIterations: [";
		unsigned int ilevel;
		for( ilevel = 0; ilevel < m_NumberOfLevels - 1; ilevel++ )
		{
			os << m_NumberOfIterations[ilevel] << ", ";
		}
		os << m_NumberOfIterations[ilevel] << "]" << std::endl;

		os << indent << "RegistrationFilter: ";
		os << m_RegistrationFilter.GetPointer() << std::endl;
		os << indent << "MovingImagePyramid: ";
		os << m_MovingImagePyramid.GetPointer() << std::endl;
		os << indent << "FixedImagePyramid: ";
		os << m_FixedImagePyramid.GetPointer() << std::endl;

		os << indent << "FieldExpander: ";
		os << m_FieldExpander.GetPointer() << std::endl;

		os << indent << "StopRegistrationFlag: ";
		os << m_StopRegistrationFlag << std::endl;

	}

	/*
	* Perform a the deformable registration using a multiresolution scheme
	* using an internal mini-pipeline
	*
	*  ref_pyramid ->  registrator  ->  field_expander --|| tempField
	* test_pyramid ->           |                              |
	*                           |                              |
	*                           --------------------------------
	*
	* A tempField image is used to break the cycle between the
	* registrator and field_expander.
	*
	*/
	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GenerateData()
	{
		// Check for NULL images and pointers
		MovingImageConstPointer movingImage = this->GetMovingImage();
		FixedImageConstPointer  fixedImage = this->GetFixedImage();
		LabelImageConstPointer  LabelImage = this->GetLabelImage();
		GradImageConstPointer MovingGradImage= this->GetMovingGradImage();
		HistogramFieldConstPointer movingHistogram=this->GetMovingHistogramField();
		HistogramFieldConstPointer fixedHistogram=this->GetFixedHistogramField();

		MultiChannelPyramidConstPointer fixedHistogramPyramid=this->GetFixedHistogramPyramid();
		MultiChannelPyramidConstPointer movingHistogramPyramid=this->GetMovingHistogramPyramid();
		// GradImageConstPointer FixedGradImage= this->GetFixedGradImage();

		if( !movingImage || !fixedImage )
		{
			itkExceptionMacro( << "Fixed and/or moving image not set" );
		}

		if( !m_MovingImagePyramid || !m_FixedImagePyramid )
		{
			itkExceptionMacro( << "Fixed and/or moving pyramid not set" );
		}

		if(!m_FixedHistogramPyramid ||!m_MovingHistogramPyramid)
		{
			itkExceptionMacro( << "Fixed and/or moving HISTOGRAM pyramid not set" );
		}

		if( !m_LabelImagePyramid || !m_LabelImagePyramid )
		{
			itkExceptionMacro( << "LabelImage pyramid not set" );
		}

		if(!m_MovingHistogramField )
		{
			itkExceptionMacro( << "Moving histogram field not set" );
		}

		if(!m_FixedHistogramField )
		{
			itkExceptionMacro( << "Fixed histogram field not set" );
		}

		if( !m_RegistrationFilter )
		{
			itkExceptionMacro( << "Registration filter not set" );
		}



		if( this->m_InitialDeformationField && this->GetInput(IDX_INIT_DEFORM) )
		{
			itkExceptionMacro( << "Only one initial deformation can be given. "
				<< "SetInitialDeformationField should not be used in "
				<< "cunjunction with SetArbitraryInitialDeformationField "
				<< "or SetInput.");
		}

		// Create the image pyramids.
		m_MovingImagePyramid->SetInput( movingImage );
		m_MovingImagePyramid->UpdateLargestPossibleRegion();

		m_FixedImagePyramid->SetInput( fixedImage );
		m_FixedImagePyramid->UpdateLargestPossibleRegion();

		
		std::cout<<"build multi-scale pyramid for moving histograms";
		m_MovingHistogramPyramid->SetBaseImage(m_MovingHistogramField);
		m_MovingHistogramPyramid->UpdateLargestPossibleRegion();
		std::cout<<"done..."<<std::endl;

		std::cout<<"build multi-scale pyramid for fixed histograms";
		m_FixedHistogramPyramid->SetBaseImage(m_FixedHistogramField);
		m_FixedHistogramPyramid->UpdateLargestPossibleRegion();
		std::cout<<"done..."<<std::endl;

		// Initializations
		m_CurrentLevel = 0;
		m_StopRegistrationFlag = false;

		unsigned int movingLevel = vnl_math_min( (int) m_CurrentLevel, 
			(int) m_MovingImagePyramid->GetNumberOfLevels() );

		unsigned int fixedLevel = vnl_math_min( (int) m_CurrentLevel, 
			(int) m_FixedImagePyramid->GetNumberOfLevels() );

		//unsigned int MovingGradLevel = vnl_math_min( (int) m_CurrentLevel, 
		//	(int) m_MovingGradImagePyramid->GetNumberOfLevels() );
		// unsigned int FixedGradLevel = vnl_math_min( (int) m_CurrentLevel, 
		//                                       (int) m_FixedGradImagePyramid->GetNumberOfLevels() );

		DeformationFieldPointer tempField = NULL;

		DeformationFieldPointer inputPtr =
			const_cast< DeformationFieldType * >( this->GetInput(0) );

		if ( this->m_InitialDeformationField )
		{
			tempField = this->m_InitialDeformationField;
		}
		else if( inputPtr )
		{
			//assert(0); //not doing this for the moment
			//let's do this
			// Arbitrary initial deformation field is set.
			// smooth it and resample

			// First smooth it
			tempField = inputPtr;

			typedef RecursiveGaussianImageFilter< DeformationFieldType,
				DeformationFieldType> GaussianFilterType;
			typename GaussianFilterType::Pointer smoother
				= GaussianFilterType::New();

			for (unsigned int dim=0; dim<DeformationFieldType::ImageDimension; ++dim)
			{
				// sigma accounts for the subsampling of the pyramid
				double sigma = 0.5 * static_cast<float>(
					m_FixedImagePyramid->GetSchedule()[fixedLevel][dim] );

				// but also for a possible discrepancy in the spacing
				sigma *= fixedImage->GetSpacing()[dim]
				/ inputPtr->GetSpacing()[dim];

				smoother->SetInput( tempField );
				smoother->SetSigma( sigma );
				smoother->SetDirection( dim );

				smoother->Update();

				tempField = smoother->GetOutput();
				tempField->DisconnectPipeline();
			}


			// Now resample
			m_FieldExpander->SetInput( tempField );

			typename FloatImageType::Pointer fi = 
				m_FixedImagePyramid->GetOutput( fixedLevel );
			m_FieldExpander->SetSize( 
				fi->GetLargestPossibleRegion().GetSize() );
			m_FieldExpander->SetOutputStartIndex(
				fi->GetLargestPossibleRegion().GetIndex() );
			m_FieldExpander->SetOutputOrigin( fi->GetOrigin() );
			m_FieldExpander->SetOutputSpacing( fi->GetSpacing());
			m_FieldExpander->SetOutputDirection( fi->GetDirection());

			m_FieldExpander->UpdateLargestPossibleRegion();
			m_FieldExpander->SetInput( NULL );
			tempField = m_FieldExpander->GetOutput();
			tempField->DisconnectPipeline();
		}

		bool lastShrinkFactorsAllOnes = false;

		while ( !this->Halt() )
		{

			if( tempField.IsNull() )
			{
				m_RegistrationFilter->SetInitialDeformationField( NULL );
			}
			else
			{
				// Resample the field to be the same size as the fixed image
				// at the current level
				m_FieldExpander->SetInput( tempField );

				typename FloatImageType::Pointer fi = 
					m_FixedImagePyramid->GetOutput( fixedLevel );
				m_FieldExpander->SetSize( 
					fi->GetLargestPossibleRegion().GetSize() );
				m_FieldExpander->SetOutputStartIndex(
					fi->GetLargestPossibleRegion().GetIndex() );
				m_FieldExpander->SetOutputOrigin( fi->GetOrigin() );
				m_FieldExpander->SetOutputSpacing( fi->GetSpacing());
				m_FieldExpander->SetOutputDirection( fi->GetDirection());

				m_FieldExpander->UpdateLargestPossibleRegion();
				m_FieldExpander->SetInput( NULL );
				tempField = m_FieldExpander->GetOutput();
				tempField->DisconnectPipeline();

				m_RegistrationFilter->SetInitialDeformationField( tempField );

			}

			// setup registration filter and pyramids 
			assert(movingLevel==fixedLevel);			//they should be the same at the moment
			m_RegistrationFilter->SetMovingImage( m_MovingImagePyramid->GetOutput(movingLevel) );
			m_RegistrationFilter->SetFixedImage( m_FixedImagePyramid->GetOutput(fixedLevel) );
			HistogramFieldType * scaledFixedHistogram=m_FixedHistogramPyramid->RetrieveLevel(movingLevel);
			m_RegistrationFilter->SetFixedHistogramField(scaledFixedHistogram);
			HistogramFieldType * scaledMovingHistogram=m_MovingHistogramPyramid->RetrieveLevel(movingLevel);
			m_RegistrationFilter->SetMovingHistogramField(scaledMovingHistogram);

			m_RegistrationFilter->SetNumberOfIterations(m_NumberOfIterations[m_CurrentLevel] );
			
			//for debuggin purpose
			//m_RegistrationFilter->SetNumberOfThreads(1);
			// cache shrink factors for computing the next expand factors.
			lastShrinkFactorsAllOnes = true;
			for( unsigned int idim = 0; idim < ImageDimension; idim++ )
			{
				if ( m_FixedImagePyramid->GetSchedule()[fixedLevel][idim] > 1 )
				{
					lastShrinkFactorsAllOnes = false;
					break;
				}
			}

			// compute new deformation field
			m_RegistrationFilter->UpdateLargestPossibleRegion();
			tempField = m_RegistrationFilter->GetOutput();
			tempField->DisconnectPipeline();

			// Increment level counter.  
			printf("level: %d\n", m_CurrentLevel);
			m_CurrentLevel++;
			movingLevel = vnl_math_min( (int) m_CurrentLevel, 
				(int) m_MovingImagePyramid->GetNumberOfLevels() );
			fixedLevel = vnl_math_min( (int) m_CurrentLevel, 
				(int) m_FixedImagePyramid->GetNumberOfLevels() );


			/*LabelLevel = vnl_math_min( (int) m_CurrentLevel, 
				(int) m_LabelImagePyramid->GetNumberOfLevels() );

			MovingGradLevel = vnl_math_min( (int) m_CurrentLevel, 
				(int) m_MovingGradImagePyramid->GetNumberOfLevels() );*/
			//   FixedGradLevel = vnl_math_min( (int) m_CurrentLevel, 
			//                                       (int) m_FixedGradImagePyramid->GetNumberOfLevels() );

			// Invoke an iteration event.
			this->InvokeEvent( IterationEvent() );

			// We can release data from pyramid which are no longer required.
			if ( movingLevel > 0 )
			{
				m_MovingImagePyramid->GetOutput( movingLevel - 1 )->ReleaseData();
				m_MovingHistogramPyramid->ReleaseLevel(movingLevel - 1);
			}
			if( fixedLevel > 0 )
			{
				m_FixedImagePyramid->GetOutput( fixedLevel - 1 )->ReleaseData();
				m_FixedHistogramPyramid->ReleaseLevel(movingLevel - 1);
			}

		} // while not Halt()

		if( !lastShrinkFactorsAllOnes )
		{
			// Some of the last shrink factors are not one
			// graft the output of the expander filter to
			// to output of this filter

			// resample the field to the same size as the fixed image
			m_FieldExpander->SetInput( tempField );
			m_FieldExpander->SetSize( 
				fixedImage->GetLargestPossibleRegion().GetSize() );
			m_FieldExpander->SetOutputStartIndex(
				fixedImage->GetLargestPossibleRegion().GetIndex() );
			m_FieldExpander->SetOutputOrigin( fixedImage->GetOrigin() );
			m_FieldExpander->SetOutputSpacing( fixedImage->GetSpacing());
			m_FieldExpander->SetOutputDirection( fixedImage->GetDirection());

			m_FieldExpander->UpdateLargestPossibleRegion();
			this->GraftOutput( m_FieldExpander->GetOutput() );
		}
		else
		{
			// all the last shrink factors are all ones
			// graft the output of registration filter to
			// to output of this filter
			this->GraftOutput( tempField );
		}

		// Release memory
		m_FieldExpander->SetInput( NULL );
		m_FieldExpander->GetOutput()->ReleaseData();
		m_RegistrationFilter->SetInput( NULL );
		m_RegistrationFilter->GetOutput()->ReleaseData();

	}


	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::StopRegistration()
	{
		m_RegistrationFilter->StopRegistration();
		m_StopRegistrationFlag = true;
	}

	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	bool
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::Halt()
	{
		// Halt the registration after the user-specified number of levels
		if (m_NumberOfLevels != 0)
		{
			this->UpdateProgress( static_cast<float>( m_CurrentLevel ) /
				static_cast<float>( m_NumberOfLevels ) );
		}

		if ( m_CurrentLevel >= m_NumberOfLevels )
		{
			return true;
		}
		if ( m_StopRegistrationFlag )
		{
			return true;
		}
		else
		{ 
			return false; 
		}

	}


	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GenerateOutputInformation()
	{

		typename DataObject::Pointer output;

		if( this->GetInput(0) )
		{
			// Initial deformation field is set.
			// Copy information from initial field.
			this->Superclass::GenerateOutputInformation();

		}
		else if( this->GetFixedImage() )
		{
			// Initial deforamtion field is not set. 
			// Copy information from the fixed image.
			for (unsigned int idx = 0; idx < 
				this->GetNumberOfOutputs(); ++idx )
			{
				output = this->GetOutput(idx);
				if (output)
				{
					output->CopyInformation(this->GetFixedImage());
				}  
			}

		}

	}


	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::GenerateInputRequestedRegion()
	{

		// call the superclass's implementation
		Superclass::GenerateInputRequestedRegion();

		// request the largest possible region for the moving image
		MovingImagePointer movingPtr = 
			const_cast< MovingImageType * >( this->GetMovingImage() );
		if( movingPtr )
		{
			movingPtr->SetRequestedRegionToLargestPossibleRegion();
		}

		// just propagate up the output requested region for
		// the fixed image and initial deformation field.
		DeformationFieldPointer inputPtr = 
			const_cast< DeformationFieldType * >( this->GetInput() );
		DeformationFieldPointer outputPtr = this->GetOutput();
		FixedImagePointer fixedPtr = 
			const_cast< FixedImageType *>( this->GetFixedImage() );

		if( inputPtr )
		{
			inputPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
		}

		if( fixedPtr )
		{
			fixedPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
		}

	}


	template <class TFixedImage, class TMovingImage, class TDeformationField, class TRealType>
	void
		LabelMultiResolutionPDEDeformableRegistration<TFixedImage,TMovingImage,TDeformationField,TRealType>
		::EnlargeOutputRequestedRegion(
		DataObject * ptr )
	{
		// call the superclass's implementation
		Superclass::EnlargeOutputRequestedRegion( ptr );

		// set the output requested region to largest possible.
		DeformationFieldType * outputPtr;
		outputPtr = dynamic_cast<DeformationFieldType*>( ptr );

		if( outputPtr )
		{
			outputPtr->SetRequestedRegionToLargestPossibleRegion();
		}

	}


} // end namespace itk

#endif
