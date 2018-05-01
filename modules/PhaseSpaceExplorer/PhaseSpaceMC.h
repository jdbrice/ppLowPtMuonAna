/*
* @Author: Daniel
* @Date:   2017-06-20 14:30:07
* @Last Modified by:   Daniel
* @Last Modified time: 2017-07-24 16:47:14
*/
#ifndef PHASE_SPACE_MC_H
#define PHASE_SPACE_MC_H


#include "HistoAnalyzer.h"
using namespace jdb;

#include "Reporter.h"
#include "vendor/loguru.h"


#include "TString.h"
#include "TRandom3.h"

#include "TLorentzVector.h"
#include "TMath.h"


class PhaseSpaceMC : public HistoAnalyzer {

protected:

	

public:
	PhaseSpaceMC() {}
	~PhaseSpaceMC() {}
	virtual void init(){
		LOG_F( INFO, "HEllo PhaseSpace" );
		
	}
	virtual void preMake(){
		book->cd();
		book->makeAll( config, nodePath + ".histograms" );
	}



	virtual void make() {
		book->cd();

		TLorentzVector lv1, lv2, lv;

		LOG_F( INFO, "Make" );
		TRandom3 rnd;
		size_t n = config.get<size_t>( "p.N", 1000 );
		for ( size_t i = 0; i < n; i++ ){

			//TOF + TOF
			lv1.SetPtEtaPhiM( rnd.Uniform( 0.17, 0.25 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv2.SetPtEtaPhiM( rnd.Uniform( 0.17, 0.25 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv = lv1+lv2;

			book->fill( "PS_TOF_TOF", lv.M(), lv.Pt() );


			//TOF + MTD
			lv1.SetPtEtaPhiM( rnd.Uniform( 0.17, 0.25 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv2.SetPtEtaPhiM( rnd.Uniform( 1.1, 11.0 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv = lv1+lv2;

			book->fill( "PS_TOF_MTD", lv.M(), lv.Pt() );


			//MTD + MTD
			lv1.SetPtEtaPhiM( rnd.Uniform( 1.1, 11.0 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv2.SetPtEtaPhiM( rnd.Uniform( 1.1, 11.0 ), rnd.Uniform( -0.5, 0.5 ), rnd.Uniform( -TMath::Pi(), TMath::Pi() ), 0.105 );
			lv = lv1+lv2;

			book->fill( "PS_MTD_MTD", lv.M(), lv.Pt() );


		}


	}


	virtual void postMake(){

		TNamed n;
		n.SetTitle( config.toXml().c_str() );
		n.Write();
	}


};



#endif