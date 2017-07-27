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

// FemtoDst
#include "FemtoDstFormat/FemtoTrackProxy.h"

// STL
#include <map>

#include "vendor/loguru.h"

class LowPtMuonFilter {
public:

	ZbRC zbUtil;
	HistoBins momentumBins;

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

	}


	map<string, XmlRange> range;

	bool pass( FemtoTrackProxy &_proxy ){

		double p = _proxy._track->mPt * cosh( _proxy._track->mEta );

		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif