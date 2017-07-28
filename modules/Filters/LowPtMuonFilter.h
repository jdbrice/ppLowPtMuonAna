/*
* @Author: Daniel
* @Date:   2017-05-05 14:22:49
* @Last Modified by:   Daniel
* @Last Modified time: 2017-05-26 10:16:12
*/
#ifndef LOW_PT_MUON_FILTER_H
#define LOW_PT_MUON_FILTER_H

#include "XmlConfig.h"
#include "XmlRange.h"
#include "HistoBins.h"
#include "XmlHistogram.h"

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// STL
#include <map>

#include "vendor/loguru.h"

class LowPtMuonFilter {
public:

	ZbRC zbUtil;
	HistoBins momentumBins;
	TFile *fitFile;
	map<string, shared_ptr<TH1> > parameters;

	LowPtMuonFilter() : zbUtil( 0.014) {}
	LowPtMuonFilter(XmlConfig &_cfg, string _nodePath) : zbUtil( 0.014) {
		load( _cfg, _nodePath );
	}
	~LowPtMuonFilter() {}

	void load( XmlConfig &_cfg, string _nodePath ){
		vector<string> rngPaths=_cfg.childrenOf( _nodePath, "Range" );

		zbUtil = _cfg.get<ZbRC>( _nodePath + ".ZbRC" );

		for ( auto p : rngPaths ){
			XmlRange rng( _cfg, p );

			range[ _cfg[ p + ":name" ] ] = rng;
			LOG_F( INFO, "[%s] = %s", _cfg[(p+":name")].c_str(), rng.toString().c_str() );
		}

		momentumBins = HistoBins(_cfg, _cfg.getString(_nodePath + ".MomentumBins") );
		LOG_F( INFO, "MomentumBins=%s", momentumBins.toString().c_str() );

		vector<string> hPaths=_cfg.childrenOf( _nodePath, "Histogram" );
		for ( auto p : hPaths ){
			XmlHistogram xhisto( _cfg, p );
			string name = _cfg.getString( p+":name" );
			parameters[ name ] =  xhisto.getTH1( );
			LOG_F( INFO, "%s : fit parameter[%s] = %p", p.c_str(), name.c_str(), parameters[ name ].get() );
		}


	}


	map<string, XmlRange> range;


	double zb( FemtoTrackProxy &_proxy, string center ){
		double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
		int pBin = momentumBins.findBin( p );
		
		double avgP     = momentumBins.bins[ pBin ] + momentumBins.binWidth( pBin ) / 2.0;
		return zbUtil.nlTof( center, _proxy._btofPid->beta(), p, avgP );
	}


	bool pass( FemtoTrackProxy &_proxy ){

		
		double p = _proxy._track->mPt * cosh( _proxy._track->mEta );
		int pBin = momentumBins.findBin( p );

		if ( pBin < 0 ) return false;
		double avgP     = momentumBins.bins[ pBin ] + momentumBins.binWidth( pBin ) / 2.0;
		double zbMu     = zbUtil.nlTof( "mu", _proxy._btofPid->beta(), p, avgP );
		double zbPi     = zbUtil.nlTof( "pi", _proxy._btofPid->beta(), p, avgP );
		
		double lambdaMu = parameters[ "mu_lambda_vs_p" ]->GetBinContent( pBin+1 );
		double lambdaPi = parameters[ "pi_lambda_vs_p" ]->GetBinContent( pBin+1 );
		double sigmaMu  = parameters[ "mu_sigma_vs_p" ]->GetBinContent( pBin+1 );
		double sigmaPi  = parameters[ "pi_sigma_vs_p" ]->GetBinContent( pBin+1 );

		zbMu = (zbMu - lambdaMu) / sigmaMu;
		zbPi = (zbPi - lambdaPi) / sigmaPi;

		if ( !range[ "nSigmaMu" ].inInclusiveRange( zbMu ) ) return false;
		if ( range[ "nSigmaPi" ].inInclusiveRange( zbPi ) ) return false;

		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif