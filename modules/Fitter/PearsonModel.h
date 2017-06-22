#ifndef TSF_PEARSON_MODEL_H
#define TSF_PEARSON_MODEL_H

// RooBarb
#include "XmlConfig.h"
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


class PearsonModel : public BaseModel
{
protected:

	

public:

	string name, yield, mu, sigma, kurt;
	double y, m, l, s;

	double e;
	string eff;

	PearsonModel( XmlConfig &cfg, string np, string _dataset = "" ) : BaseModel( cfg, np ){


		name 	= cfg.getString( np + ":name" );
		yield 	= cfg.getString( np + ":yield" );
		mu 		= cfg.getString( np + ":mu" );
		sigma 	= cfg.getString( np + ":sigma" );
		kurt	= cfg.getString( np + ":kurtosis" );

		dataset = cfg.getString( np + ":dataset", _dataset );

	}

	~PearsonModel(){}

	virtual void setVars( map< string, shared_ptr<FitVar> > &vars ){
		l = vars[ mu ]->val;
		y = vars[ yield ]->val;
		s = vars[ sigma ]->val;
		m = vars[ kurt ]->val;
	}

	virtual double eval( double x, double bw = 1.0 ){

		// e is an efficiency param
		double a1 = s * TMath::Sqrt( 2 * m - 3 );
		double a =  (y * bw) / ( a1 * TMath::Beta( m-0.5, 0.5 ) );
		double b = pow( 1 + pow( ( x - l ) / a1, 2 ), -m );

		return a * b;
	}
	double integral( double x1 = 0, double x2 = 0){

		return 0.0;
	}

	virtual string toString() { 
		return "PearsonVII( " + name + ": " + dataset + ", mu=" + dts(l) + ", sig=" +dts(s) +", y=" +dts(y) + ", kurt="+dts(m);
	}


	
};


#endif