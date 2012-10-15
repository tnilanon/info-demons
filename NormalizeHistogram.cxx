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
	std::string  outputHistogramBaseName;  /* -f option */
	std::string  histogramList;	/* -F option */
	int					nDim;						/* -d option */ //Dimension of input images, default =3

	friend std::ostream& operator<< (std::ostream& o, const arguments& args)
	{
		

		return o
			<<"Arguments structure:"<<std::endl
			<<"  Base name for the output normalized histogram field: "<<args.outputHistogramBaseName<<std::endl
			<<"  List of files of input histogram field: "<<args.histogramList<<std::endl;
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



void parseOpts (int argc, char **argv, struct arguments & args)
{
	// Command line parser
	MetaCommand command;
	command.SetParseFailureOnUnrecognizedOption( true );
	command.SetHelpCallBack(help_callback);

	// Fill some information about the software
	command.SetAuthor("Wei Liu");

	command.SetAcknowledgments("This work adapts Tom Vercauteren Demons method to distribution fields");

	command.SetDescription("Basic image registration tool with the diffeomorphic demons algorithm.");

	// Define parsing options
	command.SetOption("OutputHistogramBase","f",true,"Base name of the output histogram");
	command.SetOptionLongTag("OutputHistogramBase","fixed-image");
	command.AddOptionField("OutputHistogramBase","filename",MetaCommand::STRING,true);


	// Define parsing options for building histograms
	command.SetOption("HistogramList","F",true,"Filename containing list of histogram channels");
	command.SetOptionLongTag("HistogramList","fixed-histogram");
	command.AddOptionField("HistogramList","filename",MetaCommand::STRING,true);

	// Define parsing options for building histograms
	command.SetOption("Dimension","-d",true,"Dimensions of input images");
	command.SetOptionLongTag("Dimension","fixed-histogram");
	command.AddOptionField("Dimension","intval",MetaCommand::INT,true,"3");


	// Actually parse the command line
	if (!command.Parse(argc,argv))
	{
		exit( EXIT_FAILURE );
	}



	// Store the parsed information into a struct
	args.outputHistogramBaseName = command.GetValueAsString("OutputHistogramBase","filename");
	args.histogramList = command.GetValueAsString("HistogramList","filename");
	args.nDim=command.GetValueAsInt("Dimension","intval");
}




template <unsigned int Dimension>
void NormalizeHistograms( arguments args )
{
	// Declare the types of the images (float or double only)
	typedef float									PixelType;
	typedef itk::Image< PixelType, Dimension >		ImageType;
	typedef typename ImageType::Pointer						ImagePointer;
	typedef itk::Vector< PixelType, Dimension >		VectorPixelType;
	typedef typename itk::Image
		< VectorPixelType, Dimension >              DeformationFieldType;
	typedef itk::MultiChannelImage<ImageType>			HistogramFieldType;

	// Images we use
	typename ImageType::Pointer fixedImage = 0;
	typename ImageType::Pointer movingImage = 0;
	typename DeformationFieldType::Pointer inputDefField = 0;
	typename HistogramFieldType::Pointer pHistogram;


	// typename FixedGradImageType::Pointer fixedGradImage=0;

//	itk::HistogramPara hpara;

	pHistogram= new HistogramFieldType();
	pHistogram->LoadFromList(args.histogramList);
	std::cout<<"summing up the channels"<<std::endl;
	ImagePointer pSum=pHistogram->SumOfChannels();
	std::cout<<"normalizing...."<<std::endl;
	pHistogram->NormalizeChannels(pSum);
	std::cout<<"writing to files ..."<<std::endl;
	pHistogram->WriteToFiles(args.outputHistogramBaseName);
}


int main( int argc, char *argv[] )
{
	struct arguments args;
	parseOpts (argc, argv, args);
	std::cout<<"Starting normalizing histograms with the following arguments:"<<std::endl;
	std::cout<<args<<std::endl<<std::endl;
	std::cout<<"start running....................."<<std::endl;
	std::cout.flush();

	// Get the image dimension
	itk::ImageIOBase::Pointer imageIO;
	
	int dim=args.nDim;

	switch (dim)
	{
	case 2:
		NormalizeHistograms<2>(args);
		break;
	case 3:
		NormalizeHistograms<3>(args);
		break;
	default:
		std::cout << "Unsuported dimension" << std::endl;
		exit( EXIT_FAILURE );
	}

	return EXIT_SUCCESS;
}
