#include "HistogramField.h"

namespace itk
{

template <class MultiChannelImage > MultiChannelPyramid<MultiChannelImage>::MultiChannelPyramid()
{
	m_NumberOfLayers=2;
	m_BaseImage=0;
}

template <class MultiChannelImage > MultiChannelPyramid<MultiChannelImage>::MultiChannelPyramid(int nLevels)
{
	m_NumberOfLayers=nLevels;
	m_BaseImage=0;
}

}//namespace itk