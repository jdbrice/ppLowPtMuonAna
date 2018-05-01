#ifndef TSF_HISTO_MODEL_H
#define TSF_HISTO_MODEL_H

// RooBarb
#include "XmlConfig.h"
#include "XmlHistogram.h"
#include "Logger.h"
#include "Util.h"
using namespace jdb;

// Local
#include "Fitter/FitVar.h"
#include "Fitter/BaseModel.h"

// STL
#include <vector>
#include <map>
#include <memory>
#include <utility>

// ROOT
#include "TMath.h"
#include "TH1.h"



class HistoModel : public BaseModel
{
protected:

	

public:

	string name, yield, dataset;
	double y;

    shared_ptr<TH1> h = nullptr;
	


	PearsonModel( XmlConfig &cfg, string np, string _dataset = "" ) : BaseModel( cfg, np ){


		name 	= cfg.getString( np + ":name" );
		yield 	= cfg.getString( np + ":yield" );

        h = XmlHistogram( cfg, np ).getTH1();

		dataset = cfg.getString( np + ":dataset", _dataset );

	}

	~PearsonModel(){}

	virtual void setVars( map< string, shared_ptr<FitVar> > &vars ){
		y = vars[ yield ]->val;
	}

	virtual double eval( double x, double bw = 1.0 ){
        if (nullptr == h) return 0;
        int iBin = h->GetXaxis()->FindBin( x );
        double v = h->GetBinContent( iBin );
        double hbw = h->GetBinWidth( iBin );

        return v * hbw / bw;
	}
	double integral( double x1 = 0, double x2 = 0){
		return 0.0;
	}

	virtual string toString() { 
        return "HistoModel";
		// return "PearsonVII( " + name + ": " + dataset + ", mu=" + dts(l) + ", sig=" +dts(s) +", y=" +dts(y) + ", kurt="+dts(m);
	}


	
};


#endif