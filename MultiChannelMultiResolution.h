#pragma once
#include "HistogramField.h"
#include "itkMultiResolutionPyramidImageFilter.h"

namespace itk
{

template <class MultiChannelImage > class MultiChannelPyramid:public itk::DataObject
{
public:
	typedef MultiChannelPyramid								Self;
	typedef itk::DataObject									Superclass;
	typedef itk::SmartPointer<Self>							Pointer;
	typedef itk::SmartPointer<const Self>					ConstPointer;
	typedef typename MultiChannelImage::Pointer						MultiChannelImagePointer;
	//itkSetObjectMacro(BaseImage,MultiChannelImage);
	itkGetObjectMacro(BaseImage,MultiChannelImage);


	// use the functions provided by MultiResolutionPyramidImageFilter
	typedef typename MultiChannelImage::ImageType			InternalImageType;
	typedef typename InternalImageType::Pointer				InternalImagePointer;
	typedef MultiResolutionPyramidImageFilter<InternalImageType,InternalImageType>	SingleChannelPyramidType;
	typedef typename SingleChannelPyramidType::Pointer								SingleChannelPyramidPointer;						
	typedef std::deque<SingleChannelPyramidPointer>									ChannelPyramidQue;

	void SetBaseImage(MultiChannelImage * pt)
	{
		assert(!m_BaseImage);//currently, we require that SetBaseImage is only called once for the pyramid.
		if(pt!=m_BaseImage)
		{
			m_BaseImage=pt;
			this->Modified();
		}
		
		//now for each channel, create a multi-resolution pyramid

		int nChannels=m_BaseImage->GetNumberOfChannels();
		for(int i=0;i<nChannels;i++)
		{
			SingleChannelPyramidPointer cpp=SingleChannelPyramidType::New();
			const InternalImageType * pChannel=m_BaseImage->GetNthChannel(i); 
			cpp->SetInput(pChannel);
			cpp->SetNumberOfLevels(m_NumberOfLayers);
			_pyramids.push_back(cpp);
		}
	}
	void UpdateLargestPossibleRegion()
	{
		if(!m_BaseImage)
		{
			itkExceptionMacro(<<"Base image not set");
		}
		for(size_t i=0;i<_pyramids.size();i++)
		{
			_pyramids[i]->UpdateLargestPossibleRegion();
		}
	}
	
	MultiChannelImagePointer RetrieveLevel(int iLevel)
	{
		MultiChannelImagePointer sp=new MultiChannelImage();
		for(size_t i=0;i<_pyramids.size();i++)
		{
			SingleChannelPyramidType * cpt=_pyramids[i];
			InternalImageType * pc=cpt->GetOutput(iLevel);
			sp->PendAChannel(pc);
		}
		assert(sp->GetNumberOfChannels()==m_BaseImage->GetNumberOfChannels());
		return sp;
	}

	void ReleaseLevel(int iLevel)
	{
		for(size_t i=0;i<_pyramids.size();i++)
		{
			SingleChannelPyramidType * cpt=_pyramids[i];
			InternalImageType * pc=cpt->GetOutput(iLevel);
			pc->ReleaseData();
		}
	}

	void SetNumberOfLevels(int iLevel)
	{
		if(m_NumberOfLayers!=iLevel)
		{
			m_NumberOfLayers=iLevel;
			this->Modified();
		}
	}
	MultiChannelPyramid();
	MultiChannelPyramid(int nLevels);
protected:
	MultiChannelImagePointer		m_BaseImage;
	int								m_NumberOfLayers;
	ChannelPyramidQue				_pyramids;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MultiChannelMultiResolution.txx"
#endif