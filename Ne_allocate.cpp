#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include "fdtd3d.h"

void Ne_allocate(double *Nh, double *Re)
{
    double trash[100];

    double height(0.0);
    int bufsize = 256;
    char buf[bufsize];  //空読み用//
    int s = 25;    //読み込み初めの行//
    
    std::ifstream ifs("./IRI-2016_for_cpp/fort.7");

    if(!ifs){
        std::cout << "error" << std::endl;
        exit(0);
    }

    //空読み//
    for(int i = 0; i < s; i++){
        ifs.getline(buf, bufsize);
    }

    for(int i = 0; i < ion_L + 1; i++){

        ifs >> height >> Nh[i];
        if(Nh[i] == -1) Nh[i] = 0.0;
        else Nh[i] = Nh[i]*1.0e6;

        for(int m = 0; m < 3; m++) ifs >> trash[m];
        ifs >> Re[i];
        Re[i] = Re[i]/300.0;
        
        for(int m = 0; m < 24; m++) ifs >> trash[m];
    }

    ifs.close();
    
}