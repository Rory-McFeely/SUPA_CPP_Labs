//Header file for all my useful functions :)
#include <vector>
#include <string>

#ifndef Functions
#define Functions

//Function prototype to extract x and y data points from an input file:
std::vector<std::vector<double>> extract(std::string);


//Creating a vector to determine the magnitudes of the points in the data set:
std::vector<double> magnitude(std::vector<std::vector<double>>);


//Function which returns the coeffecients of a line of best fit using least squares method and chi squared value:
std::vector<std::vector<double>> fit_function(std::vector<std::vector<double>>, std::vector<std::vector<double>>);

//Function which returns the chi squared value of a fit:
double chi_squared(std::vector<std::vector<double>>, std::vector<std::vector<double>>, double, double);

//Function which returns x^y for each value in the data set:
std::vector<double> x_pow_y(std::vector<std::vector<double>>);

//Printing function prototypes:
void print(std::vector<std::vector<double>>);
void print(std::vector<double>, std::string);
void print(std::vector<std::vector<double>>, std::string);

//Function which takes in a user input as a selection:
int print_save_selection();

//Functions to save data to files:
void save(std::vector<std::vector<double>>, std::string);  //For saving function/chi squared
void save(std::vector<double>, std::string, std::string);  //For saving magnitude/x^y

#endif
