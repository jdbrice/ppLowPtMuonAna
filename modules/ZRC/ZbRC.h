#ifndef ZB_RECENTERING_H
#define ZB_RECENTERING_H

#include "ZRC/TofGenerator.h"

// STL
#include <map>
#include <memory>
using namespace std;

#include "XmlConfig.h"



class ZbRC {

protected:
	/**
	 * The naive sigma values for calculating the likelihood
	 */
	double tofSigma;

	shared_ptr<TofGenerator> tofGen;

	vector<string> species;

public:
	ZbRC(double tofSigma );
	~ZbRC();

	shared_ptr<TofGenerator> tofGenerator();

	// Masses
	double eMass;		// [GeV]
	double muMass;		// [GeV]
	double piMass;		// [GeV]
	
	 /**
	 * Get the expected mass from the string identifier
	 * @param  	pType 	- string identifier for species
	 * @return	mass of the given particle species
	 */
	double mass( string pType );
	
	vector<string> otherSpecies( string center );
	vector<string> allSpecies();
	vector<double> centeredTofMeans( string center, double p, vector<string> others );
	vector<double> centeredTofMeans( string center, double p );
	map<string, double> centeredTofMap( string center, double p );

	/**
	 * Calculates the difference from the expected value of 1/beta
	 * @param  	centerSpecies - particle type
	 * @param  	tof 	- the beta measured
	 * @param	p 	- momentum
	 * @return	recentered around centerSpecies
	 */
	double rTof( string centerSpecies, double beta, double p  );

	/**
	 * Calculates the difference from the expected value of 1/beta
	 * Uses the nonlinear recentering scheme
	 * @param  pType particle type
	 * @param  iHit  hit index in the nTuple
	 * @return       n Sigma from expectation
	 */
	double nlTof( string centerSpecies, double beta, double p, double avgP );

	/* Likelihood function
	 * 
	 * A gauss around the expected value with expected sigma
	 * @x     measured value
	 * @mu    expected mean
	 * @sigma expected sigma
	 * @return       the unnormalized likelihood
	 */
	static double lh( double x, double mu, double sigma );
	
};

template <>
ZbRC XmlConfig::get<ZbRC>( string path ) const;





#endif