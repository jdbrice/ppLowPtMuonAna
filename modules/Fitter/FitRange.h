#ifndef TSF_FIT_RANGE
#define TSF_FIT_RANGE

#include "Utils.h"
using namespace jdb;


class FitRange
{
public:

	string dataset;
	double min, max;

	// for dynamically tracking a variable in the fit
	string centerOn, widthFrom;
	double roiLeft, roiRight;


	FitRange( string ds, double _min, double _max, string centerOn = "", string widthFrom = "", double roi = 1.0 ){

		dataset = ds;
		min = _min;
		max = _max;

		this->centerOn = centerOn;
		this->widthFrom = widthFrom;
		this->roiLeft = roi;
		this->roiRight = roi;
	}
	FitRange( string ds, double _min, double _max, string centerOn = "", string widthFrom = "", double roiLeft = 1.0, double roiRight = 1.0 ){

		dataset = ds;
		min = _min;
		max = _max;

		this->centerOn = centerOn;
		this->widthFrom = widthFrom;
		this->roiLeft = roiLeft;
		this->roiRight = roiRight;
	}

	FitRange( const FitRange &other ){
		this->dataset = other.dataset;
		this->min = other.min;
		this->max = other.max;

		this->centerOn = other.centerOn;
		this->widthFrom = other.widthFrom;
		this->roiLeft = other.roiLeft;
		this->roiRight = other.roiRight;
	}

	
	~FitRange(){}

	bool inRange( double x ){
		if ( x < min || x > max )
			return false;
		return true;
	}

	string toString(){
		if ( "" == centerOn || "" == widthFrom )
			return "( " + dts(min) + " < " + dataset + " < " + dts(max) + " )";
		else 
			return "( " + dts(min) + " < " + dataset + " < " + dts(max) + " ) Tracking " + centerOn + " +/- " + widthFrom + " * (" + dts(roiLeft) + ", " + dts(roiRight) + ")" ;
	}
	
};	



#endif