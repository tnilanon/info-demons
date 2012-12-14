#pragma once

#include "itkImage.h"
#include "itkImageSource.h"
#include "itkWarpImageFilter.h"
#include "itkCentralDifferenceImageFunction.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkAddImageFilter.h>
#include "itkBinaryFunctorImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include <itkAccumulateImageFilter.h>
#include <deque>



namespace itk
{
	namespace Functor
	{
		template< class TInput1, class TInput2 = TInput1, class TOutput = TInput1 >
		class DivideHack
		{
		public:
			DivideHack() {}
			~DivideHack() {}
			bool operator!=(const DivideHack &) const
			{
				return false;
			}

			bool operator==(const DivideHack & other) const
			{
				return !( *this != other );
			}

			inline TOutput operator()(const TInput1 & A, const TInput2 & B) const
			{
				TOutput rs;
				if(B<1e-6)
					rs=0;
				else
					rs=A/B;
				return static_cast< TOutput >( rs);
			}
		};
	}
}

namespace itk{

	/** Type of available image forces */
	typedef enum enGradientType {
		Symmetric=0,
		Fixed=1,
		WarpedMoving=2,
		MappedMoving=3
	}GradientType;

	
template <class InternalImageType> class MultiChannelImage:public itk::DataObject
{
public:
	typedef MultiChannelImage<InternalImageType>        Self;
	typedef typename InternalImageType::Pointer				ImagePointer;
	typedef itk::DataObject									Superclass;
	typedef itk::SmartPointer<Self>							Pointer;
	typedef itk::SmartPointer<const Self>					ConstPointer;
	typedef typename InternalImageType::RegionType			RegionType;
	typedef typename InternalImageType::PixelType			PixelType;
	typedef std::deque<ImagePointer>						ImageListType;
	typedef InternalImageType								ImageType;
	itkStaticConstMacro(ImageDimension,int, InternalImageType::ImageDimension);	
	typedef typename InternalImageType::PointValueType		PointValueType;
	typedef typename InternalImageType::PointType			PointType;
	typedef	typename InternalImageType::SpacingType			SpacingType;
	typedef typename InternalImageType::DirectionType		DirectionType;

	/** Index typedef support. An index is used to access pixel values. */
	typedef typename InternalImageType::IndexType           IndexType;
	typedef typename InternalImageType::IndexValueType      IndexValueType;
	
	/** Offset typedef support. An offset represent relative position
   * between indices. */
	typedef typename InternalImageType::OffsetType          OffsetType;
	typedef typename InternalImageType::OffsetValueType		OffsetValueType;

	/** Size typedef support. A size is used to define region bounds. */
	typedef typename InternalImageType::SizeType            SizeType;
	typedef typename InternalImageType::SizeValueType       SizeValueType;

	typedef typename InternalImageType::PixelContainer		PixelContainer;
	typedef typename InternalImageType::InternalPixelType	InternalPixelType;

	///////////////////////// for warpping
	/** Interpolator type. */
	typedef double                                     CoordRepType;
	typedef InterpolateImageFunction<InternalImageType,CoordRepType> 
                                                     InterpolatorType;
	typedef typename InterpolatorType::Pointer         InterpolatorPointer;
	typedef LinearInterpolateImageFunction<InternalImageType,CoordRepType>
                                                     DefaultInterpolatorType;
	typedef typename DefaultInterpolatorType::Pointer DefaultInterpolatorPointer;
	typedef std::deque<InterpolatorPointer>				InterpolatorQue;
	/** Wapper type. */
	typedef itk::Vector< PixelType, ImageDimension >    VectorPixelType;
	typedef typename itk::Image< VectorPixelType, ImageDimension >              DeformationFieldType;
	typedef itk::WarpImageFilter<InternalImageType, InternalImageType, DeformationFieldType>           MovingImageWarpperType;
	typedef typename MovingImageWarpperType::Pointer					WarpperPointerType;
	typedef std::deque<WarpperPointerType>						WarpperQue;

	

	///////////////////////// for gradient calculating
	/** Covariant vector type. */
	typedef itk::CovariantVector<double,ImageDimension> CovariantVectorType;
	
	/** Fixed image gradient calculator type. */
	typedef itk::CentralDifferenceImageFunction<InternalImageType> GradientCalculatorType;
	typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;
	typedef std::deque<GradientCalculatorPointer>			GradientCalculatorQue;
	
	
	///////////////////////// for IO
	typedef std::deque<std::string>						ListStr;
	typedef itk::ImageFileReader< InternalImageType >	ReaderType;
	typedef typename ReaderType::Pointer				ReaderPointerType;
	typedef itk::ImageFileWriter< InternalImageType >	WriterType;
	typedef typename WriterType::Pointer				WriterPointerType;

	///////////////////////// for histogram normalizing
	typedef itk::AddImageFilter<ImageType>			AdderType;
	typedef typename AdderType::Pointer					AdderPointer;
	typedef itk::BinaryFunctorImageFilter< ImageType, ImageType, ImageType,
		itk::Functor::DivideHack<PixelType> > HackDividerType;
	typedef typename HackDividerType::Pointer		HackDividerPointer;

	typedef itk::MultiplyImageFilter<ImageType,ImageType,ImageType> MultiplierType;
	typedef typename MultiplierType::Pointer						MultiplierPointer;
	typedef itk::AccumulateImageFilter<ImageType,ImageType>			AccumulaterType;
	typedef typename AccumulaterType::Pointer						AccumulaterPointer;

	MultiChannelImage();//explicit constructor to be called on.

	void PendAChannel(ImagePointer pt)
	{
		assert(!m_IsMoving);		//not allowed to pend more channels when the histogram field starts to move.
		_channels.push_back(pt);
		SpacingType sp=pt->GetSpacing();
		if(_channels.size()>1)
		{
			assert(sp==_sp);
		}
		else
			_sp=sp;
	}
	int GetNumberOfChannels()
	{
		return _channels.size();
	}
	
	InternalImageType * GetNthChannel(int idx)
	{
		int nChannels=_channels.size();
		assert(idx>=0 && idx<nChannels);
		return _channels[idx];
	}

	/*const InternalImageType * operator [] (size_t idx)
	{
		assert(idx>=0 && idx<_channels.size());
		return _channels[idx];
	}*/

	const InternalImageType * GetFirstChannel()
	{
		return _channels[0];
	}
	// initialize warppers 
	// not needed for the static histogram field
	void InitializeMovingImages()
	{
		m_IsMoving=true;
		assert(_channels.size()>0);//
		if(m_Initialized)
			return;
		for(size_t i=0;i<_channels.size();i++)
		{
			WarpperPointerType pt=MovingImageWarpperType::New();
			DefaultInterpolatorPointer ipt=DefaultInterpolatorType::New();
			InterpolatorPointer interp=static_cast<InterpolatorType*>(ipt.GetPointer());
			_interpolators.push_back(interp);
			pt->SetInterpolator(interp);
			pt->SetEdgePaddingValue(NumericTraits<PixelType>::max());
			_warppers.push_back(pt);
		}
		m_Initialized=true;
	}
	

	// intialize gradient calculators
	// not needed for a moving histogram field, whose gradients will be calculated by hand.
	void InitializeStaticImages()
	{
		m_IsMoving=false;
		assert(_channels.size()>0);
		if(m_Initialized)
			return;
		for(size_t i=0;i<_channels.size();i++)
		{
			GradientCalculatorPointer pt=GradientCalculatorType::New();
			pt->UseImageDirectionOff();
			pt->SetInputImage(_channels[i]);
			_gradientors.push_back(pt);
		}
		m_Initialized=true;
	}

	void WarpHistogramField(const DeformationFieldType * pt)
	{
		assert(_channels.size()>0);
		assert(_channels.size()==_warppers.size()); //we create a warpper for each channel.
		ImagePointer first=_channels[0];
		PointType origin=first->GetOrigin();
		SpacingType spacing = first->GetSpacing();
		DirectionType direction = first->GetDirection();
		WarpperPointerType wpt;
		for(size_t i=0;i<_warppers.size();i++)// do a heck of warpping
		{
			wpt=_warppers[i];
			wpt->SetOutputOrigin(origin);
			wpt->SetOutputSpacing(spacing);
			wpt->SetOutputDirection(direction);
			wpt->SetInput(_channels[i]);
			wpt->SetDeformationField(pt);
			wpt->GetOutput()->SetRequestedRegion(pt->GetRequestedRegion());
			wpt->Update();
			InterpolatorPointer interp=_interpolators[i];
			interp->SetInputImage(_channels[i]);
		}
	}

	void LoadFromList(std::string fn, bool verbose=true)
	{
		ListStr ls;
		this->ReadFileList(fn, ls);
		LoadFromList(ls,verbose);
	}

	void LoadFromList(ListStr & ls, bool verbose=true)
	{
		for(size_t i=0;i<ls.size();i++)
		{
			std::string fn=ls[i];
			if(verbose)
				std::cout<<"loading file:"<<fn<<"..."<<std::endl;
			ReaderPointerType rd= ReaderType::New();
			rd->SetFileName(fn.c_str());
			rd->Update();
			ImagePointer img_ptr=rd->GetOutput();
			img_ptr->DisconnectPipeline();
			PendAChannel(img_ptr);
		}
	}

	void WriteToFiles(std::string base_fn)
	{
		char temp[255];
		std::string appendix;
		for(size_t i=0;i<_channels.size();i++)
		{
			sprintf(temp,"%04zd.hdr",i);
			appendix=temp;
			std::string fn=base_fn+appendix;
			WriterPointerType wr=WriterType::New();
			wr->SetFileName(fn.c_str());
			wr->SetInput(_channels[i]);
			wr->SetUseCompression( true );
			wr->Update();
		}
	}

	ImagePointer SumOfChannels()
	{
		ImagePointer pSum=0;
		AdderPointer pAdder=AdderType::New();
		for(size_t i=0;i<_channels.size();i++)
		{
			
			IndexType FirstIndex = _channels[i]->GetLargestPossibleRegion().GetIndex();
			SizeType ImageSize	 = _channels[i]->GetLargestPossibleRegion().GetSize();
			typename ImageType::RegionType region(FirstIndex,ImageSize);
			if(!pSum)
			{
				pSum=ImageType::New();
				pSum->SetRegions(region);
				pSum->Allocate();
				pSum->FillBuffer(0);
			}
			pAdder->SetInput1(_channels[i]);
			pAdder->SetInput2(pSum);
			pAdder->Update();
			pSum=pAdder->GetOutput();
			pAdder->GetOutput()->DisconnectPipeline();
		}
		return pSum;
	}

	void NormalizeChannels(ImagePointer pSum)
	{
		
		
		for(size_t i=0;i<_channels.size();i++)
		{
			HackDividerPointer divider=HackDividerType::New();
			ImagePointer pt=_channels[i];
			divider->SetInput1(pt);
			divider->SetInput2(pSum);
			divider->Update();
			ImagePointer dpt=divider->GetOutput();
			dpt->DisconnectPipeline();
			pt->ReleaseData();
			_channels[i]=dpt;
		}
	}

	double Coocurrance(MultiChannelImage & mc)
	{
		AccumulaterPointer acc=AccumulaterType::New();
		MultiplierPointer multiplier=MultiplierType::New();
		assert(mc._channels.size()==this->_channels.size());//
		double scooc;//sum of co-ocurrance numbers;
		for(size_t i=0;i<_channels.size();i++)
		{
			multiplier->SetInput1(_channels[i]);
			multiplier->SetInput2(mc._channels[i]);
			multiplier->Update();
			ImagePointer accRs=multiplier->GetOutput();
			for(size_t j=0;j<ImageDimension;j++){
				acc->SetInput(accRs);
				acc->SetAccumulateDirection(j);
				acc->Update();
				accRs=acc->GetOutput();
				accRs->DisconnectPipeline();
			}
			typename ImageType::IndexType sum_idx={1, 1, 1};//for 3-D images only, TOFIX
			PixelType csum=accRs->GetPixel(sum_idx);
			scooc+=csum;
		}
		return scooc;
	}

	vnl_matrix<double> CollectMovingGradient(const IndexType & idx, bool & bIsValid, vnl_vector<double> & movingValues);
	vnl_matrix<double> CollectStaticGradient(const IndexType & idx, vnl_vector<double> & fixedValues);

protected:

	bool ConsistencyCheck()//all the channels should have the same configurations. (origin, direction, spacing, size)
	{

		return true;
	}

	void ReadFileList(std::string fn, ListStr & ls)
	{
	 //read in the list of input files from the fixed group
	  std::ifstream fixedListReader;

	  fixedListReader.open(fn.c_str());

	  if(!fixedListReader.is_open())
	  {
		  std::cout << "Could not open file list:" <<fn << std::endl;
		exit( EXIT_FAILURE );
	  }

	  std::string temp;
	  fixedListReader>>temp;
	  std::cout<<"files in the list................."<<std::endl;
	  while(!fixedListReader.eof())
	  {
		ls.push_back(temp);
		std::cout<<temp<<std::endl;
		fixedListReader>>temp;
	  }
	}






	//////////////////////////////
	ImageListType	_channels;
	SpacingType		_sp;
	WarpperQue		_warppers;
	InterpolatorQue _interpolators;
	GradientCalculatorQue _gradientors; //a word made up. Don't be judgemental.
	
	bool				m_IsMoving;
	bool				m_Initialized;
};




}//namespace itk;

#ifndef ITK_MANUAL_INSTANTIATION
#include "HistogramField.txx"
#endif
