#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"
#include "geocoordinate.h"

//GEOCOORDINATE class//
void geocoordinate::set_lati(double lati){
    Lati = lati;
}

void geocoordinate::set_longi(double longi){
    Longi = longi;
}

void geocoordinate::set_alt(double alt){
    Alt = alt;
}

void geocoordinate::set_point(double lati, double longi, double alt){
    set_lati(lati);
    set_longi(longi);
    set_alt(alt);
}

void geocoordinate::geo_ijk(double lati, double longi, double alti){
    //double unit_lati = (2.0 * M_PI * R0 / 1000.0) * (lati - 135.0)/360.0;    /* 緯度1°の距離 (えびの起点) */

    //I = ;
    //J = ;
    Geo_K = alti/1.0e3;
}

void geocoordinate::set_obs(int i, int j, int k){
    Obs_I = i;
    Obs_J = j;
    Obs_K = k;
}