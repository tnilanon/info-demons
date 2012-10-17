#pragma once

#include "HistogramField.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"


namespace itk
{
	typedef struct _HistogramPara
	{
		float wSigma;			//variance parameter for soft-histogram
		int		nBins;			//number of bins
		/*range of the value in the histogram*/
		float rMax;				//maximum value
		float rMin;				//minimum value
		bool	verbose;		//print debug info?
	}HistogramPara;

	template<class TFixedImage> class HistogramBuilder
	{
	public:
		typedef MultiChannelImage<TFixedImage>							HistogramFieldType;
		typedef typename HistogramFieldType::Pointer				HistogramPointer;
		typedef std::deque<std::string>											ListStr;
		typedef TFixedImage																	ImageType;
		typedef typename ImageType::Pointer									ImagePointer;
		typedef typename TFixedImage::PixelType							PixelType;
		typedef typename TFixedImage::IndexType							IndexType;

		//for image IO
		typedef itk::ImageFileReader< ImageType >					ReaderType;
		typedef typename ReaderType::Pointer							ReaderPointer;
		typedef itk::ImageFileWriter< ImageType >					WriterType;
		typedef typename WriterType::Pointer							WriterPointer;
		
		static HistogramPointer BuildHistogramFromFileList(std::string fn, HistogramPara & para)
		{
			ListStr ls;
			ReadFileList(fn,ls);
			return BuildHistogramFromFileList(ls,para);
		}

		static bool WithinRange(IndexType curr,IndexType last)
		{
			for(int i=0;i<ImageType::ImageDimension;i++)
			{
				if(curr[i]>=last[i])
					return false;
			}
			return true;
		}
		static HistogramPointer BuildHistogramFromFileList(ListStr & f, HistogramPara & para)
		{
			//create a monster multi-channel image representing the histogram field
			HistogramPointer pHis=new HistogramFieldType();
			for(size_t i=0;i<para.nBins;i++)
			{
				ImagePointer pbin=ImageType::New();
				pHis->PendAChannel(pbin);
				//pbin->FillBuffer(0);				//Initialize to 0, important!
			}

			

			//constant preparation
			float delta_r=(para.rMax-para.rMin)/para.nBins;
			float half_delta=delta_r/2;




			assert(delta_r>1e-6);//delta_r cannot get too small;
			//
			float cutoff=4*para.wSigma;
			float wSigma2=para.wSigma*para.wSigma;
			//cutoff range of the Gaussian function is 4sigma, 
			//beyond this range, the Gaussian function value becomes ignorable.

			//load the images one by one and build histograms
			//going to be slow
			for(size_t i=0;i<f.size();i++)
			{
				std::string fn=f[i];		

				if(para.verbose)
				{
					std::cout<<"reading file:"<<fn<<std::endl;
				}

				ReaderPointer rd= ReaderType::New();
				rd->SetFileName(fn.c_str());
				rd->Update();
				ImagePointer img_ptr=rd->GetOutput();
				img_ptr->DisconnectPipeline();

				//after reading in the first image, we may allocate the memory 
				//for the histogram field
				
				IndexType FirstIndex = img_ptr->GetLargestPossibleRegion().GetIndex();
				IndexType LastIndex	= img_ptr->GetLargestPossibleRegion().GetIndex() + 
					img_ptr->GetLargestPossibleRegion().GetSize();

				typename ImageType::RegionType region(FirstIndex,img_ptr->GetLargestPossibleRegion().GetSize());
				if(!i)
				{
					for(size_t j=0;j<para.nBins;j++)
					{
						ImagePointer pbin=pHis->GetNthChannel(j);
						pbin->SetRegions(region);
						pbin->Allocate();
						pbin->FillBuffer(0);
					}
					if(para.verbose)
					{
						std::cout<<"histogram field created and initialized to 0s"<<std::endl;
					}
				}
				//typedef ImageRegionConstIterator<InputImageType> ConstIterator;
				//ConstIterator iter( img_ptr, img_ptr->GetBufferedRegion() );
				//region iterator is not suitable, as we want to establish correspondence 
				//of image pixels into the histogram field

				double sum = 0.0;
				long int count = 0;

				PixelType minValue,maxValue,meanValue;

				IndexType currIndex=FirstIndex;
				minValue = static_cast<PixelType>( img_ptr->GetPixel(currIndex) );
				maxValue = minValue;

				while(WithinRange(currIndex,LastIndex))
				{
					const PixelType value = static_cast<PixelType>( img_ptr->GetPixel(currIndex) );
					sum += static_cast<double>(value);

					if ( value < minValue ) { minValue = value; }
					if ( value > maxValue ) { maxValue = value; }

					++count;


					//now add it to the histogram.
					PixelType lower_bound=value-cutoff;
					PixelType upper_bound=value+cutoff;
					int lower_idx=(int)ceil((lower_bound-half_delta)/delta_r);
					int upper_idx=(int)floor((upper_bound-half_delta)/delta_r);
					if(lower_idx<0)
						lower_idx=0;
					assert(lower_idx<para.nBins);
					if(upper_idx>para.nBins-1)
						upper_idx=para.nBins-1;
					assert(upper_idx>=0);

					//scanned the bins within the cutoff range
					PixelType dstart=half_delta+lower_idx*delta_r-value;

					//distance between center of the starting bin to current value.
					for(int hidx=lower_idx;hidx<=upper_idx;hidx++)
					{
						PixelType weight=exp(-(dstart*dstart)/wSigma2);
						ImagePointer tpt=pHis->GetNthChannel(hidx);
						PixelType cweight=static_cast<PixelType>(tpt->GetPixel(currIndex));
						tpt->SetPixel(currIndex,cweight+weight);

						//next bin
						dstart+=delta_r;
					}



					//now increment the index
					bool fCarry=false;//check to see if there is carry
					bool fOverflow=false;
					int d=ImageType::ImageDimension-1;
					currIndex[d]++;
					if(currIndex[d]==LastIndex[d])
						fCarry=true;
					
					assert(currIndex[d]<=LastIndex[d]);
					//it is impossible for currIndex to be bigger than LastIndex

					while(fCarry)
					{
						currIndex[d]=FirstIndex[d];//set it to the smallest and carry on
						d--;
						if(d<0)//overflowed
						{
							fOverflow=true;
							break;
						}
						currIndex[d]++;
						assert(currIndex[d]<=LastIndex[d]);
						if(currIndex[d]==LastIndex[d])
							fCarry=true;
						else
							fCarry=false;
					}

					if(fOverflow)//index overflowed, last pixel has be visited.
					{
						if(para.verbose)
						{
							std::cout<<"index overflowed, current index"<<currIndex<<std::endl;
						}
						break;
					}
				}//while(currIndex<=LastIndex)

				meanValue = ( sum / static_cast<double>(count) );
				if(para.verbose)
				{
					std::cout<<"mean:"<<meanValue<<", max:"<<maxValue<<", min:"<<minValue<<std::endl;
					std::cout<<"total pixel numbers:"<<count<<std::endl;
				}

				img_ptr->ReleaseData();
			}//for(size_t i=0;i<f.size();i++)
			return pHis;
		}
	protected:
		static void ReadFileList(std::string fn, ListStr & ls)
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
	};
}
