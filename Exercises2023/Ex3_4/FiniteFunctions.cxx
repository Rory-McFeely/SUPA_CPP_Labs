#include <iostream>
#include <string>
#include <vector>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <math.h>
#include <random>

#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  //I'm guessing here RMin and RMax are being set to +/- 5.0 as default values?
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
  std::cout << "FiniteFunction constructor called..." << std::endl;
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};

/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  //ToDo write an integrator
  //Gotta get funcky wid it: Trapezoidal integration
  std::cout << "Rmin = " << m_RMin << " , Rmax = " << m_RMax << std::endl;

  double width = (m_RMax - m_RMin)/Ndiv;
  double integral = 0;
  double x0 = m_RMin;
  double area;

  for(int i=0; i<Ndiv; i++){
    double y_0 = callFunction(x0);
    double y_1 = callFunction(x0 + width);

    //Positive slope case
    if(y_0 < y_1){
      area = y_0*width + ((y_1-y_0)*width)/2;
    }
    //Negtive slope case
    else if(y_0 > y_1){
      area = y_1*width + ((y_0-y_1)*width)/2;
    }
    //Slope = 0 (edge case)
    else{
      area  = y_0*width;
    }
    //std::cout << area << std::endl;

    integral += area;
    x0 = x0 + width;
  }

  //std::cout << "The computationally determine integral = " << integral << std::endl;

  //integral = 7;
  return integral;  
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}


/*
###################
//Metropolis algorithm
###################
*/

std::vector<double> FiniteFunction::sampleFunction(int n_points){
  //Setting up the various parameters needed for random number generation:
  unsigned int seed = 100; //Setting seed
  std::mt19937 mtEngine{seed};
  std::uniform_real_distribution<double> uniformPDF{m_RMin, m_RMax};
  std::uniform_real_distribution<double> uniformPDF_2{0, 1}; //Idk if this step is necessary

  //Creating 2D vector to fill with data:
  std::vector<std::vector<double>> data;
  std::vector<double> x_vals;
  std::vector<double> y_vals;

  //First need to generate a random x point within the set range:
  double rand_x = uniformPDF(mtEngine);
  double variance = 1.0;  //Setting an arbitrary sigma

  for(int i=0; i<n_points; i++){
    //Then generate a y value using a standard distribution centred on rand_x:
    double rand_y = (1/(sqrt(variance)*sqrt(2*3.14)))*exp(-0.5*pow(((uniformPDF(mtEngine) - rand_x)/variance), 2));
    
    //std::cout << "x,y = " << rand_x << " , " << rand_y << std::endl;
    std::cout << rand_y << std::endl;

    //Next need to generate f_xi:
    double f_xi = callFunction(rand_x);

    //Performing min val test:
    double A = rand_y/f_xi;
    if(A > 1){
      A = 1;
    }

    //Generating random number T between 0 & 1:
    double T = uniformPDF_2(mtEngine);

    //Checking if T < A:
    if(T<A){
      y_vals.push_back(rand_y); //Accept y
      x_vals.push_back(rand_x);
      rand_x = rand_y;
  }

  return y_vals;
}


/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}





/* #### Gaussian Distribution #### */

//Empty constructor
NormalDistribution::NormalDistribution(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
NormalDistribution::NormalDistribution(std::vector<double> data, double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  m_variance = NULL;
  m_mean = NULL;
  m_data = data;
  this->checkPath(outfile); //Use provided string to name output files
  std::cout << "Normal Distribution constructor called..." << std::endl;
}


/*
###################
//Function eval
###################
*/ 

//Calculate mean and variance
double NormalDistribution::mean(std::vector<double> data){
  double sum = 0;
  for(int i=0; i<data.size(); i++){
    sum += data[i];
  }
  double mean = sum/data.size();
  std::cout << "Gauss mean = " << mean << std::endl;  //Print line for debugging
  m_mean = mean;
  return mean;
};

double NormalDistribution::variance(std::vector<double> data){
  double sum = 0;
  for(int i=0; i<data.size(); i++){
    sum += pow((data[i]-m_mean), 2);
  }
  double variance = sqrt(sum/data.size());
  std::cout << "Gauss variance = " << variance << std::endl;  //Print line for debugging
  m_variance = variance;
  return variance;
};

double NormalDistribution::standard_distribution(double x) {
  //Doing checks to determine if variance and mean have been calculated:
  if(m_variance == NULL){
    std::cout << "Variance not determined, calculating now" << std::endl;
    variance(m_data);
  }
  if(m_mean == NULL){
    std::cout << "Mean not detemined, calculating now" << std::endl;
    mean(m_data);
  }
  m_variance = 3.0;
  double f_x = (1/(sqrt(m_variance)*sqrt(2*3.14)))*exp(-0.5*pow(((x - m_mean)/m_variance), 2));
  return f_x;
};
double NormalDistribution::callFunction(double x) {return this->standard_distribution(x);}; //(overridable)








/* #### Cauchy-Lorentz Distribution #### */

//Empty constructor
CauchyLorentz::CauchyLorentz(){
  //Setting default values
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_gamma = 1.0;
  m_x0 = -1.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
CauchyLorentz::CauchyLorentz(std::vector<double> data, double range_min, double range_max, double mean, double gamma, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  m_data = data;
  m_gamma = gamma;
  m_x0 = mean;
  this->checkPath(outfile); //Use provided string to name output files
  std::cout << "Cauchy-Lorentz distribution constructor called..." << std::endl;
}


/*
###################
//Function eval
###################
*/ 

double CauchyLorentz::cauchy_lorentz_distribution(double x){
  //Initialising pi:
  double pi = 2*acos(0.0); //Apparently this is a good way to do it in C++ ?? Idk I just googled it
  double f_x = 1/(pi*m_gamma*((1+pow((x-m_x0)/m_gamma,2))));
  return f_x;
}

//Overriding call function:
double CauchyLorentz::callFunction(double x) {return this->cauchy_lorentz_distribution(x);};






/* #### Negative Crystal Ball Distribution #### */

//Default constructor:
CrystalBall::CrystalBall(){
  //Setting default values:
  m_RMin = -5.0;
  m_RMax = 5.0;
  m_mean = -1.0;
  m_variance = 3.0;
  m_n = 5.0;
  m_alpha = 0.5;
  this -> checkPath("DefaultFunction");
  m_Integral = NULL;
}

//Initialised constructor:
CrystalBall::CrystalBall(std::vector<double> data, double range_min, double range_max, double mean, double variance, double n, double alpha, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_data = data;
  m_mean = mean;
  m_variance = variance;
  m_n = n;
  m_alpha = alpha;
  m_Integral = NULL;
  this->checkPath(outfile);
  std::cout << "Negative crystal ball constructor called..." << std::endl;
}




/*
###################
//Function eval
###################
*/ 


double CrystalBall::CrystalBallDistribution(double x){
  //Initialising pi:
  double pi = 2*acos(0.0);

  //Intiialising return value:
  double f_x;

  //Iniialising contsants A,B,C,D,N:
  double A = (pow((m_n/abs(m_alpha)),m_n))*(exp(-(pow(abs(m_alpha),2)/2)));
  double B = (m_n/abs(m_alpha)) - abs(m_alpha);
  double C = (m_n/abs(m_alpha))*(1/(m_n-1))*(exp(-(pow(abs(m_alpha),2)/2)));
  double D = (sqrt(pi/2))*(1+erf(abs(m_alpha)/sqrt(2)));
  double N = 1/(m_variance*(C+D));

  //Using conditionals to evaluate the function based on alpha:
  double condition  = (x-m_mean)/m_variance;
  if(condition > -1*m_alpha){
    f_x = exp(-((pow((x-m_mean),2)) / (2*pow(m_variance,2))));
  }
  else{
    f_x = A*(pow((B - (x-m_mean)/m_variance), -1*m_n));
  }

  return f_x;
}


//Overriding call function::
double CrystalBall::callFunction(double x) {return this->CrystalBallDistribution(x);};
