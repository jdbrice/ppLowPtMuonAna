

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "MuonPairMaker/MuonPairMaker.h"
#include "PhaseSpaceExplorer/PhaseSpaceMC.h"
#include "Fitter/Fitter.h"
// #include "CorrectedSpectra/CorrectedSpectraMaker.h"

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);


	TaskFactory::registerTaskRunner<FemtoDstSkimmer>( "FemtoDstSkimmer" );
	TaskFactory::registerTaskRunner<MuonPairMaker>( "MuonPairMaker" );
	TaskFactory::registerTaskRunner<PhaseSpaceMC>( "PhaseSpaceMC" );
	TaskFactory::registerTaskRunner<Fitter>( "Fitter" );
	// TaskFactory::registerTaskRunner<CorrectedSpectraMaker>( "CorrectedSpectraMaker" );

	TaskEngine engine( argc, argv );

	return 0;
}
