/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelMultiResolutionPDEDeformableRegistration.h,v $
  Language:  C++
  Date:      $Date: 2009-01-26 21:45:51 $
  Version:   $Revision: 1.33 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef LabelMultiResolutionPDEDeformableRegistration_h
#define LabelMultiResolutionPDEDeformableRegistration_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "LabelDiffeomorphicDemonsRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkVectorResampleImageFilter.h"

#include "HistogramField.h"
#include "MultiChannelMultiResolution.h"

#include <vector>

namespace itk
{
/**
 * \class LabelMultiResolutionPDEDeformableRegistration
 * \brief Framework for performing multi-resolution PDE
 * deformable registration.
 *
 * LabelMultiResolutionPDEDeformableRegistration provides a generic framework
 * to peform multi-resolution deformable registration.
 *
 * At each resolution level a PDEDeformableRegistrationFilter is used
 * to register two images by computing the deformation field which will 
 * map a moving image onto a fixed image.
 *
 * A deformation field is represented as an image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the fixed image. The vector type must support element access via operator
 * []. It is assumed that the vector elements behave like floating point
 * scalars.
 *
 * The internal PDEDeformationRegistrationFilter can be set using
 * SetRegistrationFilter. By default a DemonsRegistrationFilter is used.
 *
 * The input fixed and moving images are set via methods SetFixedImage
 * and SetMovingImage respectively. An initial deformation field maybe set via
 * SetInitialDeformationField if is matches the characteristics of the coarsest
 * pyramid level. If no such assumption can be made (e.g. the deformation field
 * has the same characteristics as the input images), an initial deformation
 * field can still be set via SetArbitraryInitialDeformationField or
 * SetInput. The filter will then take care of mathching the coarsest level
 * characteristics. If no initial field is set a zero field is used as the
 * initial condition.
 *
 * MultiResolutionPyramidImageFilters are used to downsample the fixed
 * and moving images. A VectorExpandImageFilter is used to upsample
 * the deformation as we move from a coarse to fine solution.
 *
 * This class is templated over the fixed image type, the moving image type,
 * and the Deformation Field type.
 *
 * \warning This class assumes that the fixed, moving and deformation
 * field image types all have the same number of dimensions.
 *
 * \sa PDEDeformableRegistrationFilter
 * \sa DemonsRegistrationFilter
 * \sa MultiResolutionPyramidImageFilter
 * \sa VectorExpandImageFilter
 *
 * The current implementation of this class does not support streaming.
 *
 * \ingroup DeformableImageRegistration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField,  class  TRealType = float>
class LabelMultiResolutionPDEDeformableRegistration :
    public ImageToImageFilter <TDeformationField, TDeformationField>
{
public:
  /** Standard class typedefs */
  typedef LabelMultiResolutionPDEDeformableRegistration Self;
  typedef ImageToImageFilter<TDeformationField, TDeformationField>        Superclass;
  typedef SmartPointer<Self>                       Pointer;
  typedef SmartPointer<const Self>                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( LabelMultiResolutionPDEDeformableRegistration, 
                ImageToImageFilter );

  /** Fixed image type. */
  typedef TFixedImage                           FixedImageType;
  typedef typename FixedImageType::Pointer      FixedImagePointer;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointer;

  /** Moving image type. */
  typedef TMovingImage                           MovingImageType;
  typedef typename MovingImageType::Pointer      MovingImagePointer;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointer;

  /**Label image type. */
  typedef itk::Image< short, FixedImageType::ImageDimension > LabelImageType;
  typedef typename LabelImageType::Pointer		LabelImagePointer;
  typedef typename LabelImageType::ConstPointer	LabelImageConstPointer;

  /**Gradient image type */
  typedef itk::Image< float, FixedImageType::ImageDimension > GradImageType;
  typedef typename GradImageType::Pointer		GradImagePointer;
  typedef typename GradImageType::ConstPointer	GradImageConstPointer;

  /** Deformation field image type. */
  typedef TDeformationField                      DeformationFieldType;
  typedef typename DeformationFieldType::Pointer DeformationFieldPointer;

  /** ImageDimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      FixedImageType::ImageDimension);

  /** Internal float image type. */
  typedef Image<TRealType,itkGetStaticConstMacro(ImageDimension)> FloatImageType;

  /** The internal registration type. */
//  typedef PDEDeformableRegistrationFilter<FloatImageType, FloatImageType, DeformationFieldType >
  typedef LabelDiffeomorphicDemonsRegistrationFilter<FloatImageType, FloatImageType, DeformationFieldType >
                                             RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;

  /** The default registration type. */
  typedef LabelDiffeomorphicDemonsRegistrationFilter<
    FloatImageType, FloatImageType, DeformationFieldType > DefaultRegistrationType;

  /** The fixed multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<FixedImageType, FloatImageType >
                                                  FixedImagePyramidType;
  typedef typename FixedImagePyramidType::Pointer FixedImagePyramidPointer;

  /** The moving multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<MovingImageType, FloatImageType >
                                                   MovingImagePyramidType;
  typedef typename MovingImagePyramidType::Pointer MovingImagePyramidPointer;

  /** The Label multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<LabelImageType, LabelImageType>
												   LabelImagePyramidType;
  typedef typename LabelImagePyramidType::Pointer LabelImagePyramidPointer;

   /** The Gradient multi-resolution image pyramid type. */
  typedef MultiResolutionPyramidImageFilter<GradImageType, GradImageType>
												   GradImagePyramidType;
  typedef typename GradImagePyramidType::Pointer GradImagePyramidPointer;
   
  /** The deformation field expander type. */
  typedef VectorResampleImageFilter<DeformationFieldType, DeformationFieldType >
                                               FieldExpanderType;
  typedef typename FieldExpanderType::Pointer  FieldExpanderPointer;

  //********************************************************************added by wei
  //we hack the code to use as few templates as possible to avoid all those coding c*ap.
  typedef MultiChannelImage<TFixedImage> HistogramFieldType;
  typedef typename HistogramFieldType::Pointer HistogramFieldPointer;
  typedef typename HistogramFieldType::ConstPointer HistogramFieldConstPointer;
  typedef MultiChannelPyramid<HistogramFieldType>	MultiChannelPyramidType;
  typedef typename MultiChannelPyramidType::Pointer MultiChannelPyramidPointer;
  typedef typename MultiChannelPyramidType::ConstPointer MultiChannelPyramidConstPointer;
  //
  //virtual void SetFixedHistogramField(const HistogramFieldType * ptr);
  //virtual void SetMovingHistogramField(const HistogramFieldType * ptr);
  //const HistogramFieldType * GetFixedHistogramField(void) const;
  //const HistogramFieldType * GetMovingHistogramField(void) const;
  itkSetObjectMacro( FixedHistogramPyramid, MultiChannelPyramidType );
  itkSetObjectMacro( MovingHistogramPyramid, MultiChannelPyramidType );
  itkGetObjectMacro( FixedHistogramPyramid, MultiChannelPyramidType );
  itkGetObjectMacro( MovingHistogramPyramid, MultiChannelPyramidType );
  itkSetObjectMacro( FixedHistogramField, HistogramFieldType );
  itkSetObjectMacro( MovingHistogramField, HistogramFieldType );
  itkGetObjectMacro( FixedHistogramField, HistogramFieldType );
  itkGetObjectMacro( MovingHistogramField, HistogramFieldType );
  //*************wei modification ends here


  /** Set the fixed image. */
  virtual void SetFixedImage( const FixedImageType * ptr );

  /** Get the fixed image. */
  const FixedImageType * GetFixedImage(void) const;

  /** Set the moving image. */
  virtual void SetMovingImage( const MovingImageType * ptr );

  /** Get the moving image. */
  const MovingImageType * GetMovingImage(void) const;

    /** Set the label image. */
  void SetLabelImage( const LabelImageType * ptr );

  /** Get the label image. */
  const LabelImageType * GetLabelImage(void) const;

   /** Set the moving gradient image. */
  void SetMovingGradImage( const GradImageType * ptr );
  
  /** Get the moving  gradient image. */
  const GradImageType * GetMovingGradImage(void) const;
  
   /** Set the fixed gradient image. */
 // void SetFixedGradImage( const GradImageType * ptr );
  
  /** Get the fixed gradient image. */
  //const GradImageType * GetFixedGradImage(void) const;

  /** Set initial deformation field to be used as is (no smoothing, no
   *  subsampling at the coarsest level of the pyramid. */
  virtual void SetInitialDeformationField( DeformationFieldType * ptr )
    {
    this->m_InitialDeformationField=ptr;
    }

  /** Set initial deformation field. No assumption is made on the
   *  input. It will therefore be smoothed and resampled to match the
   *  images characteristics at the coarsest level of the pyramid. */
  virtual void SetArbitraryInitialDeformationField( DeformationFieldType * ptr )
    {
    this->SetInput( ptr ); 
    }

  /** Get output deformation field. */
  const DeformationFieldType * GetDeformationField(void)
  { return this->GetOutput(); }

  /** Get the number of valid inputs.  For
   * LabelMultiResolutionPDEDeformableRegistration, this checks whether the
   * fixed and moving images have been set. While
   * LabelMultiResolutionPDEDeformableRegistration can take a third input
   * as an initial deformation field, this input is not a required input.
   */
  virtual std::vector<SmartPointer<DataObject> >::size_type GetNumberOfValidRequiredInputs() const;

  /** Set the internal registrator. */
  itkSetObjectMacro( RegistrationFilter, RegistrationType );

  /** Get the internal registrator. */
  itkGetObjectMacro( RegistrationFilter, RegistrationType );
  
  /** Set the fixed image pyramid. */
  itkSetObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Get the fixed image pyramid. */
  itkGetObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Get the moving image pyramid. */
  itkGetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Set number of multi-resolution levels. */
  virtual void SetNumberOfLevels( unsigned int num );

  /** Get number of multi-resolution levels. */
  itkGetConstReferenceMacro( NumberOfLevels, unsigned int );

  /** Get the current resolution level being processed. */
  itkGetConstReferenceMacro( CurrentLevel, unsigned int );

  /** Set number of iterations per multi-resolution levels. */
  itkSetVectorMacro( NumberOfIterations, unsigned int, m_NumberOfLevels );

  /** Set the moving image pyramid. */
  itkSetObjectMacro( FieldExpander, FieldExpanderType );

  /** Get the moving image pyramid. */
  itkGetObjectMacro( FieldExpander, FieldExpanderType );

  /** Get number of iterations per multi-resolution levels. */
  virtual const unsigned int * GetNumberOfIterations() const
  { return &(m_NumberOfIterations[0]); }

  /** Stop the registration after the current iteration. */
  virtual void StopRegistration();

protected:
  LabelMultiResolutionPDEDeformableRegistration();
  ~LabelMultiResolutionPDEDeformableRegistration() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Generate output data by performing the registration
   * at each resolution level. */
  virtual void GenerateData();

  /** The current implementation of this class does not support
   * streaming. As such it requires the largest possible region
   * for the moving, fixed and input deformation field. */
  virtual void GenerateInputRequestedRegion();

  /** By default, the output deformation field has the same
   * spacing, origin and LargestPossibleRegion as the input/initial
   * deformation field.
   *
   * If the initial deformation field is not set, the output
   * information is copied from the fixed image. */
  virtual void GenerateOutputInformation();

  /** The current implementation of this class does not supprot
   * streaming. As such it produces the output for the largest
   * possible region. */
  virtual void EnlargeOutputRequestedRegion( DataObject *ptr );

  /** This method returns true to indicate that the registration should
   * terminate at the current resolution level. */
  virtual bool Halt();

private:
  LabelMultiResolutionPDEDeformableRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  RegistrationPointer        m_RegistrationFilter;
  FixedImagePyramidPointer   m_FixedImagePyramid;
  MovingImagePyramidPointer  m_MovingImagePyramid;
  LabelImagePyramidPointer   m_LabelImagePyramid;

  GradImagePyramidPointer	 m_MovingGradImagePyramid;
 // GradImagePyramidPointer	 m_FixedGradImagePyramid;
  
  FieldExpanderPointer       m_FieldExpander;
  DeformationFieldPointer    m_InitialDeformationField;

  unsigned int               m_NumberOfLevels;
  unsigned int               m_CurrentLevel;
  std::vector<unsigned int>  m_NumberOfIterations;

  /** Flag to indicate user stop registration request. */
  bool                      m_StopRegistrationFlag;

  //********************************************************************added by wei
  typename HistogramFieldType::Pointer		 m_FixedHistogramField, m_MovingHistogramField;
  MultiChannelPyramidPointer				 m_FixedHistogramPyramid, m_MovingHistogramPyramid;
  enum{IDX_INIT_DEFORM=0, IDX_FIXED_IMAGE, IDX_MOVING_IMAGE, IDX_LABEL_IMAGE, IDX_GRADE_IMAGE, 
	  IDX_FIXED_HISTOGRAM, IDX_MOVING_HISTOGRAM};
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "LabelMultiResolutionPDEDeformableRegistration.txx"
#endif


#endif
