#ifndef TSF_GAUSS_MODEL_H
#define TSF_GAUSS_MODEL_H

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


class GaussianModel : public BaseModel
{
protected:

	

public:

	string name, yield, mu, sigma;
	double y, m, s;

	double e;
	string eff;

	GaussianModel( XmlConfig &cfg, string np, string _dataset = "" ) : BaseModel( cfg, np ){


		name 	= cfg.getString( np + ":name" );
		yield 	= cfg.getString( np + ":yield" );
		mu 		= cfg.getString( np + ":mu" );
		sigma 	= cfg.getString( np + ":sigma" );
		eff 	= cfg.getString( np + ":eff", "na" );

		dataset = cfg.getString( np + ":dataset", _dataset );

	}

	GaussianModel( const GaussianModel &other ){
		this->dataset = other.dataset;
		this->name = other.name;
		this->yield = other.yield;
		this->mu = other.mu;
		this->eff = other.eff;
		this->y = other.y;
		this->m = other.m;
		this->s = other.s;
		this->s = other.s;
		this->sigma = other.sigma;
	}

	~GaussianModel(){}

	virtual void setVars( map< string, shared_ptr<FitVar> > &vars ){
		// LOG_F( INFO, "GaussianModel" );
		m = vars[ mu ]->val;
		y = vars[ yield ]->val;
		s = vars[ sigma ]->val;

		if ( vars.count( eff ) )
			e = vars[ eff ]->val; 
		else 
			e = 1.0;
	}

	virtual double eval( double x, double bw = 1.0 ){

		// e is an efficiency param
		double a = e * (y * bw) / ( s * TMath::Sqrt( 2 * TMath::Pi() ) );
		double b = TMath::Exp(  -( x - m ) * ( x - m ) / ( 2 * s * s ) );

		return a * b;
	}
	double integral( double x1 = 0, double x2 = 0){

		return 0.0;
	}

	virtual string toString() { 
		if ( "na" == eff )
			return "Gauss( " + name + ": " + dataset + ", mu=" + dts(m) + ", sig=" +dts(s) +", y=" +dts(y);
		else 
			return "Gauss( " + name + ": " + dataset + ", mu=" + dts(m) + ", sig=" +dts(s) +", y=" +dts(y) + ", e=" + dts(e);
	}


	
};


#endif