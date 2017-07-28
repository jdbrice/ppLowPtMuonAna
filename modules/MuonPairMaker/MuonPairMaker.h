#ifndef MUON_PAIR_MAKER_H
#define MUON_PAIR_MAKER_H

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
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

// Analyzers
// #include "MuonPairMaker/PairHistogramMaker.h"
// #include "MuonPairMaker/TrackHistogramMaker.h"
// #include "MuonPairMaker/MtdHistogramMaker.h"

#include "Filters/TrackFilter.h"
#include "Filters/LowPtMuonFilter.h"

#include <map>

#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"

class MuonPairMaker : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoBTofPidTraits> _rBTofPid;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	TrackFilter _trackFilter;
	LowPtMuonFilter _lowFilter;
	map<int, int> cMap;

public:
	virtual const char* classname() const {return "MuonPairMaker";}
	MuonPairMaker() {}
	~MuonPairMaker() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rHelices.setup( chain, "Helices" );
		_rBTofPid.setup( chain, "BTofPidTraits" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		_trackFilter.load( config, nodePath + ".TrackFilter" );

		cMap = config.getIntMap( nodePath + ".CentralityMap" );
		for ( int i = 0; i < 16; i++ )
			LOG_F( INFO, "cMap[%d] = %d", i, cMap[i] );

		_trackFilter.load( config, nodePath + ".TrackFilter" );
		_lowFilter.load( config, nodePath + ".LowPtMuonFilter" );
		book->cd();
	}


protected:

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();

	}
	virtual void analyzeEvent(){
	
		_event = _rEvent.get();
		book->fill( "Events", 1 );

		
		if ( cMap.count( _event->mBin16 ) == 0 ) return;
		int mappedCen = cMap[ _event->mBin16 ];
		if (mappedCen < 0) return;
		// LOG_F( INFO, "cMap[ %d ] = %d", _event->mBin16, mappedCen );
		book->fill( "mBin16", _event->mBin16 );
		book->fill( "mMappedCen", mappedCen );
		
		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;

		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rBTofPid );

			if ( _trackFilter.fail( _proxy ) ) continue;
			if ( nullptr == _proxy._btofPid ) continue;

			double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
			book->fill( "zb_p", p, _lowFilter.zb( _proxy, "mu" ) );
			if ( _lowFilter.fail( _proxy ) ) continue;

			book->fill( "zb_p_signal", p, _lowFilter.zb( _proxy, "mu" ) );

			
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
