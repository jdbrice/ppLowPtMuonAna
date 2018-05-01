/*
* @Author: Daniel
* @Date:   2017-06-20 14:30:07
* @Last Modified by:   Daniel
* @Last Modified time: 2017-07-24 16:47:14
*/

#include "HistoAnalyzer.h"
#include "HistoBins.h"
using namespace jdb;

#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
#include "Fitter/FitSchema.h"
#include "Fitter/TMinuitFitter.h"

#include "Reporter.h"
#include "RooBarb/XmlCanvas.h"
#include "TLine.h"

#include "vendor/loguru.h"


#include "TString.h"
#include "TFrame.h"

#include "RooPlotLib.h"

class Fitter : public HistoAnalyzer {

protected:

	

public:
	Fitter() {}
	~Fitter() {}
	virtual void init(){
		LOG_F( INFO, "HEllo FITTER" );
	}




	void reportFit( int pBin, TMinuitFitter &tmf, TH1 * hproj, XmlCanvas &xcanvas ){
		shared_ptr<FitSchema> _schema = tmf.schema;
		
		_schema->reportModels();

		RooPlotLib rpl;
		shared_ptr<XmlPad> xpad = xcanvas.activatePad( "fit" );
		// xpad->getRootPad()->Draw();
		LOG_F( INFO, "PAD Activated : gPad=%p", gPad );
		hproj->Draw( "he");
		hproj->GetXaxis()->SetRangeUser( -0.1, 0.2 );
		hproj->SetMinimum( 0.11 );
		hproj->SetMaximum( 1e6 );

		hproj->GetYaxis()->SetTitleOffset( 1.2 );


		gPad->Draw();
		gPad->SetLogy(0);
		gPad->Update();
		gPad->SetLogy(1);
		// TFrame *f = gPad->GetFrame();
		
		// LOG_F( INFO, "BOX(  Y: %f --> %f )", f->GetY1(), f->GetY2() );
		

		TGraph* g   = tmf.plotResult( "zb" );
		TGraph* gPi = tmf.plotResult( "zb_gPi" );
		TGraph* gMu = tmf.plotResult( "zb_gMu" );
		
		// SUM of models
		if (g) g->Draw("same");
		

		if (gPi) {
			gPi->SetLineColor(kRed);
			gPi->Draw("same");
		}
		if (gMu) {
			gMu->SetLineColor(kBlue);
			gMu->Draw("same");
		}
		
		for ( FitRange range : _schema->getRanges() ){
			if ( range.dataset != "zb" )
				continue;

			TFrame * f = gPad->GetFrame();
			LOG_F( INFO, "BOX( X: %f --> %f, Y: %f --> %f )", range.min, range.max, f->GetY1(), f->GetY2() );

			TBox * b1 = new TBox( range.min, f->GetY1(), range.max, f->GetY2() );
			b1->SetFillColorAlpha( kBlack, 0.25 );
			b1->SetFillStyle( 1001 );
			b1->Draw(  );
		}

		
		xpad = xcanvas.activatePad( "ratio" );
		

		TH1 * hratio = (TH1*)hproj->Clone( TString::Format("ratio_%d", pBin) );
		hratio->SetTitle( ";z_{b}; data / fit" );
		hratio->Reset();
		for ( int i = 1; i < hproj->GetNbinsX() + 1; i++ ){
			float x = hproj->GetBinCenter( i );
			float y = hproj->GetBinContent( i );
			

			float fv = tmf.modelEval( "zb", x );
			if ( fv < 1e-5 ) continue;
			hratio->SetBinContent( i, y / fv );
			hratio->SetBinError( i, hproj->GetBinError( i ) / y );
		}
		
		
		hratio->Draw();
		rpl.style( hratio ).set( config, "style.ratio" );
		// rpl.style( hratio ).set( "xts", 0.1 ).set( "yts", 0.1 ).set("yto", 0.5 ).set( "xto", 0.8 ).set( "yticks" ).draw();
		LOG_F( INFO, "xmin = %f", -0.1 );
		TLine *t = new TLine( -0.1, 1, 0.2, 1 );
		t->Draw("same");
		// hratio->GetYaxis()->SetRangeUser( 0.5, 1.5 );
		for ( FitRange range : _schema->getRanges() ){
			if ( range.dataset != "zb" )
				continue;
			TBox * b1 = new TBox( range.min, 0.5, range.max, 1.5 );
			b1->SetFillColorAlpha( kBlack, 0.25 );
			b1->SetFillStyle( 1001 );
			b1->Draw(  );
		}
		


		// save parameters
		for ( auto vars : _schema->getVars() ){
			TH1* h = book->get( vars.first + "_vs_p" );
			if ( nullptr == h ) continue;
			LOG_F( INFO, "[%s] = %f", vars.first.c_str(),  vars.second->val );
			h->SetBinContent( pBin, vars.second->val );
			h->SetBinError( pBin, vars.second->error );
		}
	}




	virtual void make() {

		book->cd();


		XmlCanvas xcanvas( config, "Canvas" );
		Reporter *r = new Reporter( "rp.pdf", xcanvas.getCanvas() );
		gStyle->SetOptStat(0);

		HistoBins hbP, hbEta, hbCen;
		hbP.load( config, "bins.mP" );
		hbEta.load( config, "bins.mEta" );
		hbCen.load( config, "bins.mCen" );

		

		shared_ptr<FitSchema> _schema = shared_ptr<FitSchema>( new FitSchema( config, nodePath + ".Schema" ) );
		TMinuitFitter tmf;

		tmf.schema = _schema;
		TMinuitFitter::self = &tmf;
		tmf.setupFit();

		
		TH2 * h2 = getH2D( "zb_p" );
		TH1 * hvsp = (TH1*)h2->ProjectionX("hvsp");
		
		// Create Histograms for each fit variable
		for ( auto vars : _schema->getVars() ){
			TH1D * hvarvsp = (TH1D*)hvsp->Clone(  TString(vars.first) + "_vs_p" );
			book->add( vars.first + "_vs_p", hvarvsp );
			hvarvsp->Reset();
		}

		// Setup initial conditions for fit
		tmf.fix( "mu_m" );

		
		// Loop on momentum bins
		for ( int i = 2; i < hbP.nBins(); i++ ){
			string name = "h1_" + ts( i );
			TH1 * hproj = h2->ProjectionY( name.c_str(), i+1, i+1 );

			double pmin = h2->GetXaxis()->GetBinLowEdge( i+1 );
			double pmax = pmin+h2->GetXaxis()->GetBinWidth( i+1 );
			
			LOG_F( INFO, "fitting pT bin %d = (%f, %f)", i, pmin, pmax );

			string title = dts(pmin) + "< p < " + dts(pmax) + " (GeV/c); z_{b}; dN/dz_{b}";
			hproj->SetTitle( title.c_str() );


			for ( auto vars : _schema->getVars() ){
				if ( config.exists( "Fix." + vars.first ) && i >= config.get<int>( "Fix." + vars.first ) ){
					tmf.fix( "mu_sigma" );
				}
			}

			// load dataset into fitter schema
			_schema->clearDatasets();
			_schema->loadDataset( "zb", hproj );

			// perform the fit
			tmf.fit();
			// update fitting conditions and fit again
			
			_schema->clearRanges();
			_schema->addRange( "zb", -0.2, 0.3, "pi_lambda", "pi_sigma", config.get<float>( "pi:roiLeft", 3 ), config.get<float>( "pi:roiRight", 3 ) );
			_schema->addRange( "zb", -0.2, 0.3, "mu_lambda", "mu_sigma", config.get<float>( "mu:roiLeft", 1 ), config.get<float>( "mu:roiRight", 1 ) );
			
			_schema->updateRanges( );
			
			// final fit
			tmf.fit();
			
			reportFit( i, tmf, hproj, xcanvas );
			r->savePage();
			r->saveImage( "plot_fit_and_ratio_" + ts(i) + ".pdf" );
			if ( 5 == i ){
				r->saveImage( "plots.C" );
			}
			
		}

		xcanvas.cd();
		RooPlotLib rpl;
		for ( auto vars : _schema->getVars() ){
			TH1* h = book->get( vars.first + "_vs_p" );
			if ( nullptr == h ) continue;
			h->SetTitle( TString(  "; p (GeV/c);"+vars.first ) );
			h->Draw();
			rpl.style( h ).set( config, "style.fit_param" ).draw();
			r->savePage();
		}


		r->close();
		delete r;


	}


	virtual void postMake(){

		TNamed n;
		n.SetTitle( config.toXml().c_str() );
		n.Write();
	}


};