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
    int nBins = 49;

    //For Cauchy-Lorentz distribution  {2.6 seems to be good value}
    double gamma = 2.6;  

    //Crystal Ball variables: (These actually give a pretty good fit, but obv not as good as Gaussian)
    double mean = -1.0;
    double variance = 3.0;
    double n = 5.0;
    double alpha = 2.0;

    FiniteFunction finite_distribution = FiniteFunction(range_min, range_max, "Finite_function");
    NormalDistribution normal_distribution = NormalDistribution(data, range_min, range_max, "Normal_distribution");
    CauchyLorentz cauchy_distribution = CauchyLorentz(data, range_min, range_max, mean, gamma, "Cauchy_Lorentz_distribution");
    CrystalBall crystal_ball_distribution = CrystalBall(data, range_min, range_max, mean, variance, n, alpha, "Crystal_ball_distribution");

    //Plotting:
    finite_distribution.plotData(data, nBins, true);
    finite_distribution.plotFunction();

    normal_distribution.plotData(data, nBins, true);
    normal_distribution.plotFunction();

    cauchy_distribution.plotData(data, nBins, true);
    cauchy_distribution.plotFunction();

    crystal_ball_distribution.plotData(data, nBins, true);
    crystal_ball_distribution.plotFunction();


}