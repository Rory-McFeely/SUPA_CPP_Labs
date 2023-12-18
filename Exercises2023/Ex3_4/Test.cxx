//Author: Rory McFeely
//Date: 04/12/23

#include <iostream>
#include <string>
#include "FiniteFunctions.h"
#include <fstream>


int main(){
    //Extracting the random data set
    std::string line;
    std::string fileName = "Outputs//data//MysteryData01441.txt";
    std::ifstream inputFile(fileName);

    //Creating a vector to fill with data points:
    std::vector<double> data;
    while(std::getline(inputFile, line)){
        data.push_back(std::stod(line));
        //std::cout << line << std::endl;   //Print line for debugging
    }

    //Declaring variables
    double range_min = -10.0;
    double range_max = 10.0;
    double gamma = 2.6;  //For Cauchy-Lorentz distribution  {2.6 seems to be good value}
    int nBins = 49;
    std::string output_file = "testOutput";

    //Crystal Ball variables: (These actually give a pretty good fit, but obv not as good as Gaussian)
    double mean = -1.0;
    double variance = 3.0;
    double n = 5.0;
    double alpha = 2.0;

    //FiniteFunction test = FiniteFunction(range_min, range_max, output_file);
    //NormalDistribution test = NormalDistribution(data, range_min, range_max, output_file);
    //CauchyLorentz test = CauchyLorentz(data, range_min, range_max, gamma, output_file);
    CrystalBall test = CrystalBall(data, range_min, range_max, mean, variance, n, alpha, output_file);
    test.plotData(data, nBins, true);
    test.plotFunction();


}