#include "HistogramField.h"

#include "LabelESMDemonsRegistrationFunction.h"

namespace itk{



template <class InternalImageType> MultiChannelImage<InternalImageType>::MultiChannelImage():m_Initialized(false)
{
	m_IsMoving=false;
}

template <class InternalImageType> vnl_matrix<double> 
MultiChannelImage<InternalImageType>::CollectStaticGradient(const IndexType & idx,vnl_vector<double> & staticValues)
{
	if(!ConsistencyCheck())
	{
		itkExceptionMacro(<<"Channels of the histogram field image must be consistent");
	}

	assert(!m_IsMoving); //must not be called on a moving histogram field
	size_t nChannel=_channels.size();
	assert(nChannel>0);
	vnl_matrix<double> d(nChannel,ImageDimension);
	staticValues=vnl_vector<double>(nChannel);

	//get a bunch of metrix and boundary values
	ImagePointer first=_channels[0];
	PointType origin=first->GetOrigin();
	SpacingType spacing = first->GetSpacing();
	DirectionType direction = first->GetDirection();
	IndexType FirstIndex = first->GetLargestPossibleRegion().GetIndex();
    IndexType LastIndex =  first->GetLargestPossibleRegion().GetIndex() + 
     _channels[0]->GetLargestPossibleRegion().GetSize();

	//std::cout<<"FirstIndex of static gradient"<<FirstIndex<<std::endl;//for debugging and verification purpose
	//std::cout<<"LastIndex of static gradient"<<LastIndex<<std::endl;

	for(size_t i=0;i<nChannel;i++)
	{
		const CovariantVectorType fixedGradient=_gradientors[i]->EvaluateAtIndex(idx);
		staticValues[i]=_channels[i]->GetPixel(idx);
		for(int j=0;j<ImageDimension;j++)
			d(i,j)=fixedGradient[j];
	}
	return d;
}

template <class InternalImageType> vnl_matrix<double> 
MultiChannelImage<InternalImageType>::CollectMovingGradient(const IndexType & index, bool & bIsValid, 
															vnl_vector<double> & movingValues)
{
	assert(m_IsMoving); //must be called on a moving histogram field
	size_t nChannel=_channels.size();
	assert(nChannel>0);
	vnl_matrix<double> d(nChannel,ImageDimension);
	bIsValid=true;
	movingValues=vnl_vector<double>(nChannel);

	ImagePointer first=_channels[0];
	PointType origin=first->GetOrigin();
	SpacingType spacing = first->GetSpacing();
	DirectionType direction = first->GetDirection();
    IndexType FirstIndex = _channels[0]->GetLargestPossibleRegion().GetIndex();
    IndexType LastIndex = _channels[0]->GetLargestPossibleRegion().GetIndex() + 
     _channels[0]->GetLargestPossibleRegion().GetSize();

	//std::cout<<"FirstIndex of moving gradient"<<FirstIndex<<std::endl;
	//std::cout<<"LastIndex of moving gradient" <<LastIndex<<std::endl;

	for(size_t c=0;c<nChannel;c++)				//loop over every channel and calculate the gradient vector 
	{
		
		//calculate by hand, due to some numeric nuances happened during the warping step.
		WarpperPointerType wpt=_warppers[c];


		// Get moving image related information
		// check if the point was mapped outside of the moving image using
		// the "special value" NumericTraits<MovingPixelType>::max()
		PixelType movingPixValue
			= wpt->GetOutput()->GetPixel( index );
		
 

		if( movingPixValue == NumericTraits <PixelType>::max() )
		{
			movingValues[c]=0;
			bIsValid=false;
			return d;
		}
		
		//if(bIsValid)
		//	std::cout<<"catch one"<<std::endl;

		const double movingValue = static_cast<double>( movingPixValue );
		
		movingValues[c]=movingValue;

		// We compute the gradient more or less by hand.
		// We first start by ignoring the image orientation and introduce it afterwards 
		//CovariantVectorType usedOrientFreeGradientTimes2;
  
		/*if( (m_UseGradientType==Symmetric) || 
			(m_UseGradientType==WarpedMoving) )

		{*/
			// we don't use a CentralDifferenceImageFunction here to be able to
			// check for NumericTraits<MovingPixelType>::max()
		CovariantVectorType warpedMovingGradient;
		IndexType tmpIndex = index;
		for( unsigned int dim = 0; dim < ImageDimension; dim++ )
		{
			// bounds checking
			if( FirstIndex[dim]==LastIndex[dim]
			|| index[dim]<FirstIndex[dim]
			|| index[dim]>=LastIndex[dim] )
			{
				warpedMovingGradient[dim] = 0.0;//out of bound, no gradients
				continue;
			}
			else if ( index[dim] == FirstIndex[dim] ) //at the front boundary
			{
				// compute derivative
				tmpIndex[dim] += 1;
				movingPixValue = wpt->GetOutput()->GetPixel( tmpIndex );
				if( movingPixValue == NumericTraits <PixelType>::max() )
				{
					// weird crunched border case
					warpedMovingGradient[dim] = 0.0;
				}
				else
				{
					// forward difference
					warpedMovingGradient[dim] = static_cast<double>( movingPixValue ) - movingValue;
					warpedMovingGradient[dim] /= spacing[dim]; 
				}
				tmpIndex[dim] -= 1;
				continue;
			}
			else if ( index[dim] == (LastIndex[dim]-1) )
			{
				// compute derivative
				tmpIndex[dim] -= 1;
				movingPixValue = wpt->GetOutput()->GetPixel( tmpIndex );
				if( movingPixValue == NumericTraits<PixelType>::max() )
				{
					// weird crunched border case
					warpedMovingGradient[dim] = 0.0;
				}
				else
				{
					// backward difference
					warpedMovingGradient[dim] = movingValue - static_cast<double>( movingPixValue );
					warpedMovingGradient[dim] /= spacing[dim]; 
				}
				tmpIndex[dim] += 1;
				continue;
			}//else if ( index[dim] == (LastIndex[dim]-1) )//at the rear boundary


			// current pixel is not at any boundaries
			tmpIndex[dim] += 1;
			movingPixValue = wpt->GetOutput()->GetPixel( tmpIndex );
			if ( movingPixValue == NumericTraits
				<PixelType>::max() ) //if invalid, try neighbors on the other side
			{
				// backward difference
				warpedMovingGradient[dim] = movingValue;

				tmpIndex[dim] -= 2;
				movingPixValue = wpt->GetOutput()->GetPixel( tmpIndex );
				if( movingPixValue == NumericTraits<PixelType>::max() )
				{
					// weird crunched border case
					warpedMovingGradient[dim] = 0.0; //both neighbors are bad, give up.
				}
				else
				{
					// backward difference
					warpedMovingGradient[dim] -= static_cast<double>(
						wpt->GetOutput()->GetPixel( tmpIndex ) );

					warpedMovingGradient[dim] /= spacing[dim];
				}
			}
			else // the neighbor is perfectly ok.
			{
				warpedMovingGradient[dim] = static_cast<double>( movingPixValue );

				tmpIndex[dim] -= 2;
				movingPixValue = wpt->GetOutput()->GetPixel( tmpIndex );
				//see if the other neighbor is also ok. If Yes, then use both of them to calculate the gradient.
				if ( movingPixValue == NumericTraits<PixelType>::max() )
				{
					//the second neighbor is not good, do a one-sided difference
					// forward difference
					warpedMovingGradient[dim] -= movingValue;
					warpedMovingGradient[dim] /= spacing[dim];
				}
				else
				{
					// normal case, central difference
					warpedMovingGradient[dim] -= static_cast<double>( movingPixValue );
					warpedMovingGradient[dim] *= 0.5 / spacing[dim];
				}
			}
			tmpIndex[dim] += 1;//shift tmpIndex back to normal
		}// for( unsigned int dim = 0;
		
		//now stuff the gradient vector into vnl_matrix

		for(int j=0;j<ImageDimension;j++)
		{
			d(c,j)=warpedMovingGradient[j];
		}
	}//for(size_t c=0;c<nChannel;c++) //calculating gradients for each channel

	//note that the calculated gradients here are orient-free. Orientation will be added later.
	return d;
}


}