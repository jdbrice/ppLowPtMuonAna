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
#include "Filters/MuonBDTFilter.h"
#include "Filters/MuonMLPFilter.h"

#include <map>

#include "TLorentzVector.h"

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

	TrackFilter _lowPtTrackFilter;
	TrackFilter _mtdTrackFilter;
	LowPtMuonFilter _lowFilter;
	MuonBDTFilter _bdtFilter;
	MuonMLPFilter _mlpFilter;

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

		_lowPtTrackFilter.load( config, nodePath + ".LowPtTrackFilter" );
		_mtdTrackFilter.load( config, nodePath + ".MtdTrackFilter" );
		_lowFilter.load( config, nodePath + ".LowPtMuonFilter" );
		_bdtFilter.load( config, nodePath + ".MuonBDTFilter" );
		_mlpFilter.load( config, nodePath + ".MuonMLPFilter" );
		book->cd();
	}


protected:

	vector<FemtoTrackProxy> pos_mtd;
	vector<FemtoTrackProxy> neg_mtd;
	vector<FemtoTrackProxy> pos_tof;
	vector<FemtoTrackProxy> neg_tof;

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();
	}

	virtual void makePairs( vector<FemtoTrackProxy> &col1, vector<FemtoTrackProxy> &col2, string prefix ){
		TLorentzVector lv1, lv2, lv;
		for ( FemtoTrackProxy& _proxy1 : col1 ){
			for ( FemtoTrackProxy& _proxy2 : col2 ){
				lv1.SetPtEtaPhiM( _proxy1._track->mPt, _proxy1._track->mEta, _proxy1._track->mPhi, 0.105 );
				lv2.SetPtEtaPhiM( _proxy2._track->mPt, _proxy2._track->mEta, _proxy2._track->mPhi, 0.105 );
				lv = lv1 + lv2;

				book->fill( prefix + "_pt_mass", lv.M(), lv.Pt() );
			} // loop col1
		} // loop col2
	}


	virtual void makePairs( vector<FemtoTrackProxy> &col1, string prefix ){
		TLorentzVector lv1, lv2, lv;

		for ( size_t i = 0; i < col1.size(); i++ ){
			for ( size_t j = i; j < col1.size(); j++ ){
				if ( i == j ) continue;
				FemtoTrackProxy& _proxy1 = col1[i];
				FemtoTrackProxy& _proxy2 = col1[j];
				lv1.SetPtEtaPhiM( _proxy1._track->mPt, _proxy1._track->mEta, _proxy1._track->mPhi, 0.105 );
				lv2.SetPtEtaPhiM( _proxy2._track->mPt, _proxy2._track->mEta, _proxy2._track->mPhi, 0.105 );
				lv = lv1 + lv2;

				book->fill( prefix + "_pt_mass", lv.M(), lv.Pt() );
			}
		}
	}

	virtual void analyzeEvent(){
	
		_event = _rEvent.get();
		book->fill( "Events", 1 );

		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;

		size_t nTOF = 0;
		size_t nMTD = 0;

		pos_mtd.clear();
		neg_mtd.clear();
		pos_tof.clear();
		neg_tof.clear();

		for (size_t i = 0; i < nTracks; i++ ){
			_proxy.assemble( i, _rTracks, _rHelices, _rBTofPid );
			_proxy.setMtdPidTraits( _rMtdPid );

			int charge = _proxy._track->charge();
			
			if ( _proxy._track->mBTofPidTraitsIndex >= 0 && _lowPtTrackFilter.pass( _proxy ) ){
				double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
				double zb = _lowFilter.zb( _proxy, "mu" );
				book->fill( "zb_p", p, zb );
				if ( _lowFilter.pass( _proxy ) ){
					book->fill( "zb_p_signal", p,zb );
					if ( charge > 0 )
						pos_tof.push_back( _proxy );
					else 
						neg_tof.push_back( _proxy );
					nTOF++;
				}
			}
			
			if ( _mtdTrackFilter.pass( _proxy ) &&  _proxy._track->mMtdPidTraitsIndex >= 0
				&& fabs(_proxy._mtdPid->mDeltaY) < 10
				&& fabs(_proxy._mtdPid->mDeltaZ) < 10
				&& fabs(_proxy._mtdPid->mDeltaTimeOfFlight) < 0.5
				){

				// DeltaZ > 60 is HACK until new MC is trained
				// float bdt = _bdtFilter.evaluate( _proxy );
				// float mlp = _mlpFilter.evaluate( _proxy );
				// book->fill( "BDT", bdt );
				// book->fill( "MLP", mlp );
				// book->fill( "MLP_vs_BDT", bdt, mlp );

				// book->fill( "dY_BDT", bdt, _proxy._mtdPid->mDeltaY );
				// book->fill( "dZ_BDT", bdt, _proxy._mtdPid->mDeltaZ );
				// book->fill( "dTOF_BDT", bdt, _proxy._mtdPid->mDeltaTimeOfFlight );
				// book->fill( "nSigmaPi_BDT", bdt, _proxy._track->nSigmaPion() );
				// book->fill( "nHitsFit_BDT", bdt, fabs(_proxy._track->mNHitsFit) );
				// book->fill( "dca_BDT", bdt, _proxy._track->gDCA() );
				// book->fill( "cell_BDT", bdt, _proxy._mtdPid->cell() );
				// book->fill( "module_BDT", bdt, _proxy._mtdPid->module() );
				// book->fill( "backleg_BDT", bdt, _proxy._mtdPid->backleg() );

				// book->fill( "dY_MLP", mlp, _proxy._mtdPid->mDeltaY );
				// book->fill( "dZ_MLP", mlp, _proxy._mtdPid->mDeltaZ );
				// book->fill( "dTOF_MLP", mlp, _proxy._mtdPid->mDeltaTimeOfFlight );
				// book->fill( "nSigmaPi_MLP", mlp, _proxy._track->nSigmaPion() );
				// book->fill( "nHitsFit_MLP", mlp, fabs(_proxy._track->mNHitsFit) );
				// book->fill( "dca_MLP", mlp, _proxy._track->gDCA() );
				// book->fill( "cell_MLP", mlp, _proxy._mtdPid->cell() );
				// book->fill( "module_MLP", mlp, _proxy._mtdPid->module() );
				// book->fill( "backleg_MLP", mlp, _proxy._mtdPid->backleg() );

				if ( charge > 0 )
					pos_mtd.push_back( _proxy );
				else 
					neg_mtd.push_back( _proxy );
				nMTD++;

			}
		} // loop on tracks

		book->fill( "nTof_vs_nMtd", nMTD, nTOF );
		book->fill( "nPos_vs_nNeg", pos_tof.size() + pos_mtd.size(), neg_tof.size() + neg_mtd.size() );

		makePairs( pos_tof, neg_mtd, "uls" );
		makePairs( neg_tof, pos_mtd, "uls" );
		makePairs( pos_tof, pos_mtd, "ls" );
		makePairs( neg_tof, neg_mtd, "ls" );

		makePairs( neg_tof, "tof_ls" );
		makePairs( pos_tof, "tof_ls" );
		makePairs( pos_tof, neg_tof, "tof_uls" );


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
