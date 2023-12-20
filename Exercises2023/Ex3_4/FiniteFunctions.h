#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF


//1/(1+x**2) Class:
class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  double sampleFunction(double x);  //Function implementing metropolis algorithim

  //Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  
private:
  double invxsquared(double x); //The default functional form
};




//Normal distribution class:
class NormalDistribution : public FiniteFunction{

  public:
    NormalDistribution();
    NormalDistribution(std::vector<double> data, double range_min, double range_max, std::string outfile);
    virtual double callFunction(double x); //Call the function with value x (Overridable)
    double variance(std::vector<double> data);
    double mean(std::vector<double> data);
  
  protected:
    double m_mean;
    double m_variance;
    std::vector<double> m_data;

  private:
    double standard_distribution(double x);

};





//Cauchy-Lorentz distribution class:
class CauchyLorentz : public FiniteFunction{

  public:
    CauchyLorentz();
    CauchyLorentz(std::vector<double> data, double range_min, double range_max, double mean, double gamma, std::string outfile);
    virtual double callFunction(double x); //Call the function with value x (Overridable)
  
  protected:
    double m_gamma;
    double m_x0;
    std::vector<double> m_data;

  private:
    double cauchy_lorentz_distribution(double x);

};




//Negative Crystal Ball distribution class:
class CrystalBall : public FiniteFunction{

  public:
    CrystalBall(); //Default constructor
    CrystalBall(std::vector<double> data, double range_min, double range_max, double mean, double variance, double n, double alpha, std::string outfile); //Custom constructor
    virtual double callFunction (double x);

  protected:
    double m_mean;
    double m_variance;
    double m_n;
    double m_alpha;
    std::vector<double> m_data;  //I should have just implemented this into FiniteFunction class and then i wouldn't have to redeclare it in each child class :(

  private:
    double CrystalBallDistribution(double x);  //Actual crystal ball distribution function
};