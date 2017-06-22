#ifndef TSF_POLY_MODEL_H
#define TSF_POLY_MODEL_H

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


class PolyModel : public BaseModel
{
protected:

	

public:

	string name;
	vector<string> coeff_names;
	vector<double> coeff;
	int n = 0;

	PolyModel( XmlConfig &cfg, string np, string _dataset = "" ) : BaseModel( cfg, np ){

		LOG_F( INFO, "PolyModel" );

		name 	= cfg.getString( np + ":name" );
		coeff_names 	= cfg.getStringVector( np + ":c" );
		dataset = cfg.getString( np + ":dataset", _dataset );
		n = coeff_names.size();

	}

	PolyModel( const PolyModel &other ){
		this->coeff_names.clear();
		this->coeff.clear();
		for ( int i = 0; i < other.coeff_names.size(); i++ ){
			this->coeff_names.push_back( other.coeff_names[i] );
			this->coeff.push_back( other.coeff[i] );
		}
	}

	~PolyModel(){}

	virtual void setVars( map< string, shared_ptr<FitVar> > &vars ){
		// LOG_F( INFO, "SetVars" );
		coeff.clear();
		// LOG_F( INFO, "n = %d", coeff_names.size() );

		for ( int i = 0; i < coeff_names.size(); i++ ){

			double v = vars[ coeff_names[i] ]->val;
			// LOG_F( INFO, "var[%s] = %f", coeff_names[i].c_str(), v );
			coeff.push_back( v );
		}
		// for ( string cname : coeff_names ){
		// 	LOG_F( INFO, "var[%s] = %f", cname.c_str(), vars[ cname ]->val );
		// 	coeff.push_back( vars[ cname ]->val );
		// }
		// LOG_F( INFO, "DONE" );
	}

	virtual double eval( double x, double bw = 1.0 ){

		double v = 0;
		for ( int i = 0; i < n; i++ ){
			v += coeff[i] * pow( x, i );
		}
		return v;
	}
	double integral( double x1 = 0, double x2 = 0){

		return 0.0;
	}

	virtual string toString() { 
		string m = "Polynomial( ";
		string sep="";
		for ( int i = 0; i < coeff.size(); i++ ){
			m += sep + dts( coeff[i] ) + " x^" + ts(i);
			sep = "+";
		}
		return m + ")";
	}


	
};


#endif