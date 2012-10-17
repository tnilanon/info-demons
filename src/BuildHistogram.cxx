/**
* This is a small tool that shows how to use the diffeomorphic demons algorithm.
* The user can choose if diffeomorphic, additive or compositive demons should be used.
* The user can also choose the type of demons forces, or other parameters;
*
* \author Tom Vercauteren, INRIA & Mauna Kea Technologies
*/


#include <itkCommand.h>

#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <itkFastSymmetricForcesDemonsRegistrationFilter.h>
#include <itkGridForwardWarpImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <itkTransformFileReader.h>
#include <itkTransformToDeformationFieldSource.h>
#include <itkVectorCentralDifferenceImageFunction.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include <itkWarpHarmonicEnergyCalculator.h>
#include <itkWarpImageFilter.h>
#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "LabelDiffeomorphicDemonsRegistrationFilter.h"
#include "LabelMultiResolutionPDEDeformableRegistration.h"
#include "LabelMultiResolutionPDEDeformableRegistration.h"
#include "itkGradientMagnitudeImageFilter.h"
#include <metaCommand.h>
#include "HistogramField.h"
#include "SoftHistogram.h"

#include <errno.h>
#include <iostream>
#include <limits.h>



struct arguments
{
	std::string  fixedGroupBaseName;  /* -o option */
	std::string  fixedGroupList;			/* -i option */
	float wSigma;									/* -w option */ //Window size for soft-histogram
	int nBins;										/* -b					*/ //number of bins in the histogram
	//min and max of the value
	float rMin;										/* -m option*/
	float rMax;										/* -M option*/

	friend std::ostream& operator<< (std::ostream& o, const arguments& args)
	{


		return o
			<<"Arguments structure:"<<std::endl
			<<"  Output base name for the computed histogram field: "<<args.fixedGroupBaseName<<std::endl
			<<"  Input list of images: "<<args.fixedGroupList<<std::endl
			<<"  Width of the Gaussian window function: "<<args.wSigma<<std::endl
			<<"	 Number of bins in the histogram: "<<args.nBins<<std::endl
			<<"	 Range of the histogram value: "<<args.rMin<<" to "<<args.rMax<<std::endl;

	}
};

void help_callback()
{
	std::cout<<std::endl;
	std::cout<<"Copyright (c) 2008 INRIA and Mauna Kea Technologies"<<std::endl;
	std::cout<<"Code: Tom Vercauteren"<<std::endl;
	std::cout<<"Report bugs to <tom.vercauteren \\at maunakeatech.com>"<<std::endl;

	exit( EXIT_FAILURE );
};

int atoi_check( const char * str )
{
	char *endptr;
	long val= strtol(str, &endptr, 0);

	/* Check for various possible errors */
	if ( (errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
		|| val>=INT_MAX || val<=INT_MIN )
	{
		std::cout<<std::endl;
		std::cout<<"Cannot parse integer. Out of bound."<<std::endl;
		exit( EXIT_FAILURE );
	}

	if (endptr == str || *endptr!='\0')
	{
		std::cout<<std::endl;
		std::cout<<"Cannot parse integer. Contains non-digits or is empty."<<std::endl;
		exit( EXIT_FAILURE );
	}

	return val;
}


std::vector<unsigned int> parseUIntVector( const std::string & str)
{
	std::vector<unsigned int> vect;

	std::string::size_type crosspos = str.find('x',0);

	if (crosspos == std::string::npos)
	{
		// only one uint
		vect.push_back( static_cast<unsigned int>( atoi_check(str.c_str()) ));
		return vect;
	}

	// first uint
	vect.push_back( static_cast<unsigned int>(
		atoi_check( (str.substr(0,crosspos)).c_str()  ) ));

	while(true)
	{
		std::string::size_type crossposfrom = crosspos;
		crosspos =  str.find('x',crossposfrom+1);

		if (crosspos == std::string::npos)
		{
			vect.push_back( static_cast<unsigned int>(
				atoi_check( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str()  ) ));
			return vect;
		}

		vect.push_back( static_cast<unsigned int>(
			atoi_check( (str.substr(crossposfrom+1,crosspos-crossposfrom-1)).c_str()  ) ));
	}
}


void parseOpts (int argc, char **argv, struct arguments & args)
{
	// Command line parser
	MetaCommand command;
	command.SetParseFailureOnUnrecognizedOption( true );
	command.SetHelpCallBack(help_callback);

	// Fill some information about the software
	command.SetAuthor("Wei Liu");

	command.SetAcknowledgments("doodle");

	command.SetDescription("Build soft-histograms for Jensen-Shannon Demons.");

	// Define parsing options
	command.SetOption("OutputHistogramBase","o",true,"Base name of the output histogram");
	command.SetOptionLongTag("OutputHistogramBase","histogram-basename");
	command.AddOptionField("OutputHistogramBase","filename",MetaCommand::STRING,true);

	// Define parsing options for building histograms
	command.SetOption("InputImageList","i",true,"Filename containing list of input images for building fixed histogram field");
	command.SetOptionLongTag("InputImageList","input-list");
	command.AddOptionField("InputImageList","filename",MetaCommand::STRING,true);

	command.SetOption("HistogramWindowWidth","w",false,"Variance of the Gaussian window functon");
	command.SetOptionLongTag("HistogramWindowWidth","gauss-window-width");
	command.AddOptionField("HistogramWindowWidth","floatval",MetaCommand::FLOAT,true,"4.0");

	command.SetOption("NumberOfBins","b",false,"Number of bins in the histogram");
	command.SetOptionLongTag("NumberOfBins","bin-number");
	command.AddOptionField("NumberOfBins","intval",MetaCommand::INT,true,"128");

	command.SetOption("RangeMinimum","m",false,"Minimum value of the histogram ");
	command.SetOptionLongTag("RangeMinimum","range-minimum");
	command.AddOptionField("RangeMinimum","floatval",MetaCommand::FLOAT,true,"0");

	command.SetOption("RangeMaximum","M",false,"Maximum value of the histogram");
	command.SetOptionLongTag("RangeMaximum","range-maximum");
	command.AddOptionField("RangeMaximum","floatval",MetaCommand::FLOAT,true,"255");//command.SetOption("RangeOfHistogramValues","i",false," Range of the histogram value< MinxMax>");
  //command.SetOptionLongTag("RangeOfHistogramValues","histogram-range");
  //command.AddOptionField("RangeOfHistogramValues","uintfloat",MetaCommand::STRING,true,"0x255");



	// Actually parse the command line
	if (!command.Parse(argc,argv))
	{
		exit( EXIT_FAILURE );
	}



	// Store the parsed information into a struct
	args.fixedGroupBaseName = command.GetValueAsString("OutputHistogramBase","filename");
	args.fixedGroupList = command.GetValueAsString("InputImageList","filename");
	args.wSigma =command.GetValueAsFloat("HistogramWindowWidth","floatval");
	args.rMin =command.GetValueAsFloat("RangeMinimum","floatval");
	args.rMax =command.GetValueAsFloat("RangeMaximum","floatval");
	args.nBins = command.GetValueAsInt("NumberOfBins", "intval");
}




template <unsigned int Dimension>
void BuildHistograms( arguments args )
{
	// Declare the types of the images (float or double only)
	typedef float									PixelType;
	typedef short									LabelPixelType;
	typedef itk::Image< PixelType, Dimension >		ImageType;
	typedef itk::Image< LabelPixelType, Dimension > LabelImageType;

	typedef itk::Vector< PixelType, Dimension >		VectorPixelType;
	typedef typename itk::Image
		< VectorPixelType, Dimension >              DeformationFieldType;
	typedef itk::MultiChannelImage<ImageType>			HistogramFieldType;

	// Images we use
	typename ImageType::Pointer fixedImage = 0;
	typename ImageType::Pointer movingImage = 0;
	typename LabelImageType::Pointer LabelImage =0;
	typename DeformationFieldType::Pointer inputDefField = 0;
	typename HistogramFieldType::Pointer fixedHistogramField;
//	typename HistogramFieldType::Pointer movingHistogramField;

	typedef float										GradPixelType;
	typedef  itk::Image<GradPixelType,Dimension>		MoveGradImageType;
	typename MoveGradImageType::Pointer					movingGradImage=0;
	// typename FixedGradImageType::Pointer fixedGradImage=0;

	itk::HistogramPara hpara;
	hpara.wSigma =args.wSigma;
	hpara.nBins=args.nBins;
	hpara.rMax=args.rMax;
	hpara.rMin=args.rMin;
	hpara.verbose=true;
	fixedHistogramField= itk::HistogramBuilder<ImageType>::BuildHistogramFromFileList(args.fixedGroupList,hpara);
	fixedHistogramField->WriteToFiles(args.fixedGroupBaseName);
}


int main( int argc, char *argv[] )
{
	struct arguments args;
	parseOpts (argc, argv, args);
	std::cout<<"Starting demons registration with the following arguments:"<<std::endl;
	std::cout<<args<<std::endl<<std::endl;
	std::cout<<"start running....................."<<std::endl;
	std::cout.flush();
	// FIXME uncomment for debug only
	// itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

	// Get the image dimension
	itk::ImageIOBase::Pointer imageIO;
	int dim=3;

	switch (dim)
	{
	case 2:
		BuildHistograms<2>(args);
		break;
	case 3:
		BuildHistograms<3>(args);
		break;
	default:
		std::cout << "Unsuported dimension" << std::endl;
		exit( EXIT_FAILURE );
	}

	return EXIT_SUCCESS;
}
