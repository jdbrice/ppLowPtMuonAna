/*
* @Author: Daniel
* @Date:   2017-06-20 15:15:15
* @Last Modified by:   Daniel
* @Last Modified time: 2017-06-21 16:23:47
*/

#ifndef TMINUIT_FITTER_H
#define TMINUIT_FITTER_H
#include "Fitter/FitSchema.h"

#include "TMinuit.h"
#include "TGraph.h"

#include <memory>

#include "vendor/loguru.h"

class TMinuitFitter
{
public:

	shared_ptr<FitSchema> schema;

	// ordered list of parameters
	// these are the parameter names in minuit parameter order
	vector<string> parNames;

	// the minuit interface
	unique_ptr<TMinuit> minuit;

	TMinuitFitter() {}
	~TMinuitFitter() {}

	static TMinuitFitter *self;

	static void tminuitFCN(int &npar, double *gin, double &f, double *par, int flag){
		double fnVal = 0.0;
		updateParameters( npar, par );

		bool useRange = self->schema->constrainFitRange();
		string method = self->schema->getMethod();
		// double norm = self->getNorm();

		// loop on datasets
		for ( auto k : self->schema->datasets ){

			string ds = k.first;

			// loop on datapoints
			for ( auto d : k.second ){
				
				// test for ranges
				if ( useRange && !self->schema->inRange( ds, d.x )){
					continue;
				}


				if ( "chi2" == method){
					// Minimize by chi^2
					if ( 0 == d.ey  )
						continue;
					double modelVal = self->modelEval( ds, d.x );
					fnVal += chi2( d.y, modelVal, d.ey );	
				} else if ( "poisson" == method ){
					// Minimize by poisson errors ie error = sqrt( N )
					if ( 0 >= d.y  )
						continue;
					double modelVal = self->modelEval( ds, d.x );
					fnVal += poisson( d.y, modelVal );	
				} else if ( "fractional" == method ){
					// Minimize by fractional errors ie error =  N 
					if ( 0 >= d.y || 0 >= d.ey )
						continue;
					double modelVal = self->modelEval( ds, d.x );
					fnVal += fractional( d.y, modelVal );	
				} else {
					cout << "No Fit method" << endl;
				}
				
			} // loop on data points
		}

		f = (fnVal );
	}

	void fit(){

		double arglist[10];
		arglist[ 0 ] = 50000;
		arglist[ 1 ] = 1.0;
		int iFlag = -1;
		string status = "na";

		minuit->mnexcm( "MINI", arglist, 1, iFlag );
		status = minuit->fCstatu;


		updateParameters();
	}

	void setupFit(){
		minuit = unique_ptr<TMinuit>( new TMinuit( schema->numParams() ) );
		
		for ( auto k : schema->getVars() ){
			if ( k.second->exclude ) continue;
			parNames.push_back( k.first );
		}

		for ( int i = 0; i < parNames.size(); i++ ){
			if (schema->var( parNames[ i ] )->exclude) continue;

			if ( 0 <= schema->var( parNames[ i ] )->error )
				schema->var( parNames[ i ] )->error = 0.001;

			minuit->DefineParameter( i, 						// parameter index
						parNames[ i ].c_str(), 					// name
						schema->var( parNames[ i ] )->val, 		// initial value
						schema->var( parNames[ i ] )->error,	// intial step size
						schema->var( parNames[ i ] )->min,		// limit min
						schema->var( parNames[ i ] )->max );	// limit max

			if ( schema->var( parNames[ i ] )->fixed )
				minuit->FixParameter( i );

		}

		minuit->SetPrintLevel( schema->getVerbosity() );
		minuit->SetFCN( tminuitFCN );

		// fitIsGood = false;
	}

	static double chi2( double data, double pred, double ey ){
		double q = (data - pred) / ey;
		return q*q;
	};
	static double poisson( double data, double pred){
		double q = (data - pred) * (data - pred) / data;
		return q;
	};

	static double fractional( double data, double pred){
		double q = (data - pred) * (data - pred) / (data * data);
		return q;
	};

	void nop(  );

	static double modelEval( string dataset, double x ){
		// LOG_F( INFO, "ModelEval( %s, %f )", dataset.c_str(), x );
		// cout << "s is " << self->schema << endl;
		// cout << "s is " << self->schema << endl;
		// cout << "s is " << self->schema << endl;
		
		// cout << "is: " << self->schema->datasets.count( dataset ) << endl;
		double bw = self->schema->datasets[ dataset ].pointNear( x ).bw;
		double val = 0.0;
		for ( auto k : self->schema->models ){
			// cout << "name : " << k.first << ", dataset=" << k.second->dataset << endl;
			if ( k.second->dataset == dataset ){ // second part is looking for active

				val += k.second->eval( x, bw );
			}
		}
		return val;
	}

	static void updateParameters( int npar =0, double * pars = 0){
			
		// update the variables
		for ( int i = 0; i < self->parNames.size(); i++ ){
			
			if ( npar >= self->parNames.size() ){
				double val = pars[ i ];
				self->schema->var( self->parNames[ i ] )->val = val;

			} else {
				double val = 0, err = 0;
				self->minuit->GetParameter( i, val, err );
				self->schema->var( self->parNames[ i ] )->val = val;
				self->schema->var( self->parNames[ i ] )->error = err;

			}
		}
		self->schema->updateModels( );
	}

	TGraph * plotResult( string datasetOrModel ){
		LOG_F( INFO, "Plotting : %s", datasetOrModel.c_str() );

		string datasetName = "";
		string modelName = "";

		shared_ptr<BaseModel> gm;
		bool found = true;

		auto dsIt = schema->datasets.find( datasetOrModel );
		auto mIt = schema->models.find( datasetOrModel );

		LOG_F( INFO, "finded?" );

		if ( dsIt != schema->datasets.end() ){ // show all models for this dataset
			LOG_F( INFO, "found dataset %s", datasetOrModel.c_str() );
			datasetName = dsIt->first;
		} else if ( mIt != schema->models.end() ){ // show just this model
			LOG_F( INFO, "found model %s", datasetOrModel.c_str() );
			gm = mIt->second;
			modelName = mIt->first;
			datasetName = mIt->second->dataset;
		} else {
			found = false;
		}

		if ( !found ){
			return new TGraph();
		}

		// gm->eval( 0, 0 );
		// if ( "" != modelName && !self->players[ modelName ] ){
		// 	return new TGraph();
		// }

		pair<double, double> xBounds = schema->datasets[ datasetName ].rangeX(  );

		LOG_F( INFO, "Plotting %s from ( %f, %f )", datasetName.c_str(), xBounds.first, xBounds.second );	

		double stepsize = (xBounds.second - xBounds.first) / 500.0;
		vector<double> vx, vy;

		if ( found ){

			int i = 0;
			for ( double x = xBounds.first; x < xBounds.second; x+= stepsize ){
				// LOG_F( INFO, "x=%f", x );
				
				
				FitDataPoint fdp = schema->datasets[ datasetName ].pointNear( x );
				
				
				vx.push_back(x);
				if ( !gm ){
					vy.push_back( modelEval( datasetName, x ) );
				}
				else 
					vy.push_back( gm->eval( x, fdp.bw ) );

				i++;
			}	


		}		

		TGraph * graph = new TGraph( vx.size(), vx.data(), vy.data() );
		return graph;
	}
	
};



#endif