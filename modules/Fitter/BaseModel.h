/*
* @Author: Daniel
* @Date:   2017-06-21 15:11:08
* @Last Modified by:   Daniel
* @Last Modified time: 2017-06-21 16:16:16
*/

#ifndef TSF_BASE_MODEL_H
#define TSF_BASE_MODEL_H

// RooBarb
#include "XmlConfig.h"
#include "Util.h"
using namespace jdb;

// Local
#include "Fitter/FitVar.h"

// STL
#include <vector>
#include <map>
#include <memory>
#include <utility>

// ROOT
#include "TMath.h"

#include "vendor/loguru.h"

class BaseModel
{
protected:

	

public:
	string dataset;

	BaseModel(){}

	BaseModel( XmlConfig &cfg, string np){

	}
	
	BaseModel( const BaseModel &other ){

	}

	~BaseModel(){}

	virtual void setVars( map< string, shared_ptr<FitVar> > &vars ){
		LOG_F( INFO, "BaseModel" );
	}

	virtual double eval( double x, double bw = 1.0 ){
		return 1.0;
	}
	
	double integral( double x1 = 0, double x2 = 0){
		return 0.0;
	}

	virtual string toString() { 
		return "BaseMode";
	}


	
};


#endif