#include "ZRC/ZbRC.h"


template <>
ZbRC XmlConfig::get<ZbRC>( string path ) const{
	ZbRC zb( get<float>( path + ":sigma" ) );
	return zb;
}

ZbRC::ZbRC( double tofSigma ){
	this->tofSigma = tofSigma;

	tofGen = shared_ptr<TofGenerator>(  new TofGenerator() );

	species = { "e", "mu", "pi" };

	// in GeV / c^2
	eMass 		= 0.000510998;
	muMass		= 0.1056583715; 
	piMass 		= 0.1395702;
}
ZbRC::~ZbRC(){
}

shared_ptr<TofGenerator> ZbRC::tofGenerator() { return tofGen; }

double ZbRC::mass( string pType ){
	if ( "pi" == pType )
		return piMass;
	if ( "mu" == pType )
		return muMass;
	if ( "e" == pType )
		return eMass;
	return -10.0;	
}
vector<string> ZbRC::otherSpecies( string center ){
	vector<string> res;
	for ( string plc : species ){
		if ( plc != center )
			res.push_back( plc );
	}
	return res;
}
vector<string> ZbRC::allSpecies(){
	return species;
}
vector<double> ZbRC::centeredTofMeans( string center, double p, vector<string> others ){

	double cMean = tofGen->mean( p, mass( center ) );
	
	vector<double> res;
	for ( string plc : others ){
		double m = (tofGen->mean( p, mass( plc ) ) - cMean);
		res.push_back( m );
	}

	return res;
}
vector<double> ZbRC::centeredTofMeans( string center, double p ){

	double cMean = tofGen->mean( p, mass( center ) );
	
	vector<double> res;
	for ( string plc : species ){
		double m = (tofGen->mean( p, mass( plc ) ) - cMean);
		res.push_back( m );
	}

	return res;
}
map<string, double> ZbRC::centeredTofMap( string center, double p ){

	double cMean = tofGen->mean( p, mass( center ) );
	
	map<string, double> res;
	for ( string plc : species ){
		double m = (tofGen->mean( p, mass( plc ) ) - cMean);
		res[ plc ] = m;
	}

	return res;
}

double ZbRC::rTof( string centerSpecies, double beta, double p  ){
	return ( ( 1.0 / beta ) - tofGen->mean( p, mass( centerSpecies ) ) );
}

double ZbRC::nlTof( string centerSpecies, double beta, double p, double avgP ){

	const double sigma = tofSigma; 
	const double tof = 1.0 / beta;
	
	// mean for this species
	//const double mu =  tofGen->mean( p, mass( centerSpecies ) );
	const double muAvg =  tofGen->mean( avgP, mass( centerSpecies ) );
	
	double n1 = 0;
	double d1 = 0;

	for ( string plc : species ){
		
		double iMu =  tofGen->mean( p, mass( plc ) ) ;
		double iMuAvg =  tofGen->mean( avgP, mass( plc ) ) ;
		
		
		double iL = lh( tof, iMu, sigma );
		
		double w = tof + iMuAvg - iMu;
		
		n1 += (iL * w);
		d1 += iL;

	}

	if ( 0 == n1 || 0 == d1){
		return tof - muAvg;
	}
	
	double nTof = (n1/d1) - ( muAvg );
	
	return nTof;

}

double ZbRC::lh( double x, double mu, double sigma ){

	double a = sigma * sqrt( 2 * TMath::Pi() );
	double b = ( x - mu );
	double c = 2 * sigma*sigma;
	double d = (1/a) * exp( -b*b / c );

	return d;
}