/*
* @Author: Daniel
* @Date:   2017-06-20 14:30:07
* @Last Modified by:   Daniel
* @Last Modified time: 2017-06-21 17:41:35
*/

#include "HistoAnalyzer.h"
#include "HistoBins.h"
using namespace jdb;

#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "Fitter/FitSchema.h"
#include "Fitter/TMinuitFitter.h"

#include "Reporter.h"

#include "vendor/loguru.h"




class Fitter : public HistoAnalyzer {

protected:

	

public:
	Fitter() {}
	~Fitter() {}
	virtual void init(){
		LOG_F( INFO, "HEllo FITTER" );
	}

	virtual void make() {

		HistoBins hbP, hbEta, hbCen;
		hbP.load( config, "bins.mP" );
		hbEta.load( config, "bins.mEta" );
		hbCen.load( config, "bins.mCen" );

		shared_ptr<FitSchema> _schema = shared_ptr<FitSchema>( new FitSchema( config, nodePath + ".Schema" ) );
		TMinuitFitter tmf;

		tmf.schema = _schema;
		// _schema->addRange( "zb", -0.2, 0.3 );
		// _schema->addRange( "zb", -0.2, 0.3, "pi_lambda", "pi_sigma", 1 );
		// _schema->addRange( "zb", -0.2, 0.3, "mu_lambda", "mu_sigma", 1 );
		TMinuitFitter::self = &tmf;
		tmf.setupFit();

		Reporter *r = new Reporter("rp.pdf", 1600, 1024 );

		TH2 * h2 = getH2D( "zb_p" );
		for ( int i = 2; i < hbP.nBins(); i++ ){
			LOG_F( INFO, "fitting pT bin %d", i );
			string name = "h1_" + ts( i );
			TH1 * hproj = h2->ProjectionY( name.c_str(), i+1, i+1 );

			_schema->clearDatasets();
			_schema->loadDataset( "zb", hproj );
			double pmin = h2->GetXaxis()->GetBinLowEdge( i+1 );
			double pmax = pmin+h2->GetXaxis()->GetBinWidth( i+1 );
			string title = dts(pmin) + "< p < " + dts(pmax) + " (GeV/c); #beta^{-1}; dN/d#beta^{-1}";
			r->newPage();
			hproj->Draw();
			hproj->SetTitle( title.c_str() );
			tmf.fit();
			tmf.fit();

			TGraph* g=tmf.plotResult( "zb" );
			TGraph* gPi=tmf.plotResult( "zb_gPi" );
			TGraph* gMu=tmf.plotResult( "zb_gMu" );
			
			if (g) g->Draw("same");
			gPad->SetLogy(1);
			if (gPi) {
				gPi->SetLineColor(kRed);
				gPi->Draw("same");
			}
			if (gMu) {
				gMu->SetLineColor(kBlue);
				gMu->Draw("same");
			}
			_schema->updateRanges();
			_schema->reportModels();
			
			r->savePage();
		}


		r->close();
		delete r;


	}


};