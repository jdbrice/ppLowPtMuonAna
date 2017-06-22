#ifndef FEMTO_DST_SKIMMER_H
#define FEMTO_DST_SKIMMER_H

#include "TreeAnalyzer.h"
#include "HistoBins.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"
#include "FemtoDstFormat/FemtoBTofPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

// Analyzers
// #include "FemtoDstSkimmer/PairHistogramMaker.h"
#include "FemtoDstSkimmer/TrackHistogramMaker.h"
// #include "FemtoDstSkimmer/MtdHistogramMaker.h"

#include "Filters/TrackFilter.h"

#include "FemtoDstSkimmer/TofGenerator.h"

#include <map>

#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"

class FemtoDstSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoBTofPidTraits> _rBTofPid;

	TrackHistogramMaker thm;
	bool trackQA = false;
	TrackFilter _trackFilter;

	HistoBins hbPt, hbEta, hbCen;
	map<int, int> cMap;

	TofGenerator tofGen;

	double rTof( double zb, double p, double m ){
		double zbp = zb - tofGen.mean( p, m );
		return zbp;
	}


public:
	virtual const char* classname() const {return "FemtoDstSkimmer";}
	FemtoDstSkimmer() {}
	~FemtoDstSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rHelices.setup( chain, "Helices" );
		_rBTofPid.setup( chain, "BTofPidTraits" );


		trackQA = config.getBool( "TrackHistogramMaker:enable", true );
		if ( trackQA ) {
			LOG_F( INFO, "Creating Track QA" );
			thm.setup( config, "TrackHistogramMaker", book );
		}
		_trackFilter.load( config, nodePath + ".TrackFilter" );
		book->cd();


		hbPt.load( config, "bins.signal.mP" );
		hbEta.load( config, "bins.signal.mEta" );
		hbCen.load( config, "bins.signal.mCen" );

		cMap = config.getIntMap( nodePath + ".CentralityMap" );
		for ( int i = 0; i < 16; i++ )
			LOG_F( INFO, "cMap[%d] = %d", i, cMap[i] );

		


	}

	static string charge_string( int c ){
		if ( c < 0 )
			return "n";
		if ( c > 0 )
			return "p";
		return "E";
	}
	static string bin_name( int cen, int eta, int c ){
		stringstream sstr;
		sstr << "cen_" << cen << "_eta_" << eta << "_c_" << charge_string( c );
		return sstr.str(); 
	}


protected:

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();

		book->cd();


		//for ( int iPt = 0; iPt < hbPt.nBins(); iPt++ ){
			for ( int iEta = 0; iEta < hbEta.nBins(); iEta++ ){
				for ( int iCen = 0; iCen < hbCen.nBins(); iCen++ ){
					for ( int iCharge : { -1, 1 } ) {
						string bn = bin_name( iCen, iEta, iCharge );
						LOG_F( INFO, "Making histo : %s", bn.c_str() );

						book->clone( "zd_p", "zd_p_" + bn );
						book->clone( "zb_p", "zb_p_" + bn );

						// book->clone( "yLocal_p", "yLocal_p_" + bn );
						// book->clone( "zLocal_p", "zLocal_p_" + bn );
					}
				}
			}
		// }

	}
	virtual void analyzeEvent(){
	
		_event = _rEvent.get();
		book->fill( "Events", 1 );

		book->fill( "mBin16", _event->mBin16 );

		int mappedCen = cMap[ _event->mBin16 ];
		
		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;
		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rBTofPid );

			if ( nullptr == _proxy._btofPid ) continue;

			if ( _trackFilter.fail( _proxy ) ) continue;
			if ( trackQA ) thm.analyze( _proxy );

			int iEta = hbEta.findBin( _proxy._track->mEta );
			if ( iEta < 0 ) continue;
			string tname = bin_name( mappedCen, iEta, _proxy._track->charge() );


			if ( fabs(_proxy._btofPid->yLocal()) > 1.6 ) continue;
			if ( fabs(_proxy._btofPid->zLocal()) > 2.8 ) continue;

			double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
			double zb = rTof( 1.0/_proxy._btofPid->beta(), p, 0.105658374 );

			book->fill( "zd_p", p, _proxy._track->nSigmaPion() );
			book->fill( "zb_p", p, zb );

			book->fill( "zd_p_" + tname, p, _proxy._track->nSigmaPion() );
			book->fill( "zb_p_" + tname, p, zb );

			book->fill( "yLocal_p", p, _proxy._btofPid->yLocal() );
			book->fill( "zLocal_p", p, _proxy._btofPid->zLocal() );

		}
	}


	virtual void postEventLoop(){
		TreeAnalyzer::postEventLoop();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
		
	}
	
};

#endif
