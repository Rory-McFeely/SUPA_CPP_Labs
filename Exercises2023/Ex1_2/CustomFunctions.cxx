//These functions may or may not be useful
//Rory McFeely

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include "CustomFunctions.h"



//The function below takes a filename/directory as an argument and extracts x and y data points, returning them as a 2D array
std::vector<std::vector<double>> extract(std::string FileName){

    std::string line;
    std::string fileName = FileName;
    std::ifstream inputFile(fileName);
    std::string delimeter = ",";
    int nline = 0;

    //Initialising two vectors to fill with x and y data:
    std::vector<double> x_points;
    std::vector<double> y_points;

    std::vector<std::vector<double>> data;

        //This bit checks that the file opened succesfully and that nothing got fucked along the way
    if (!inputFile.is_open()) {
        std::cout << "Error opening file: " << fileName << std::endl;
        //Note: You can't declare a vector/double etc. within the return line so you have to do it in a line above
        std::vector<std::vector<double>> empty_vec;
        return empty_vec;
    }

    //So long as file opens properly the while loop below runs through the opened file and prints each line to the terminal
    else{
        while(std::getline(inputFile, line)){
            //The if statement below checks to see if the loop is on the header of the text file and skips it
            if(nline == 0){
                nline +=1;
                continue;
            }

            //The line below is used to create a string containing the characters up to the first instance of the delimiter. In this case the x data point
            //.subst creates a substring, and .find() will find the position of the delimiter variable value in the string being analysed
            std::string x_val_temp = line.substr(0, line.find(delimeter));
            std::string y_val_temp = line.substr(line.find(delimeter)+1, -1);
            //std::cout << x_val_temp << " , " << y_val_temp << std::endl;
            x_points.push_back(std::stod(x_val_temp));
            y_points.push_back(std::stod(y_val_temp));

            
            nline+=1;
        }
        data.push_back(x_points);
        data.push_back(y_points);
        }

    return data;
}


//Function which takes a vector of vectors as an input and prints their values:
void print(std::vector<std::vector<double>> Data){
    //Checking if data contains > 5 lines:
    int nlines = Data[0].size();
    if (nlines>5){
        std::cout << "Number of lines in data set: " << nlines << "\n" << "Only 5 lines will be printed" << std::endl;
    }
    //Obviously making the assumption that data is ordered as  x,y in the array:
    //Could do some stuff with formatting here to make them print in nice even collumns
    std::cout << "  x     ,     y" << std::endl;

    for(int i=0; i<5; i++){
        std::cout << Data[0][i] << " ,  " << Data[1][i] << "\n";
    }
}

//Overloading the print function so it can also take a 1D vector i.e. magnitudes:
void print(std::vector<double> data, std::string header){
    std::cout << header << std::endl;
    for(int i=0; i < data.size(); i++){
        std::cout << data[i] << std::endl;
    }
}

//Overloading the print function to work with the fit function
void print(std::vector<std::vector<double>> data, std::string header){
    std::cout << header << std::endl;
    std::cout << "y = " << data[0][0] << "x + " << data[0][1] << std::endl;
    std::cout << "Chi squared = " << data[1][0] << std::endl;
}


//Function which takes a vector of vectors as an input and determines the magnitude of each x,y in relation to an origin at (0,0)
std::vector<double> magnitude(std::vector<std::vector<double>> Data){
    //Initialising a vector to fill with the magnitude values
    std::vector<double> magnitudes;

    for(int i=0; i<Data[0].size(); i++){
        double temp_magnitude = sqrt(pow(Data[0][i], 2) + pow(Data[1][i], 2));
        magnitudes.push_back(temp_magnitude);
    }

    return magnitudes;
}


//Fucntion which uses least squares method to determine the linear best fit coefficients:
std::vector<std::vector<double>> fit_function(std::vector<std::vector<double>> Data, std::vector<std::vector<double>> error){
    //Initialising variables to store the sums, sum of squares etc.
    double sum_x = 0;
    double sum_y = 0;
    double sum_xy = 0;
    double sum_x_square = 0;

    //Initialising the two fitting variables:
    double p;
    double q;

    //Initialising return vector:
    std::vector<std::vector<double>> best_fits;

    //Loops through each data point and adds to the sum of x,y, product of xy and x^2 variables
    for(int i=0; i<Data[0].size(); i++){
        sum_x += Data[0][i];
        sum_y += Data[1][i];
        sum_xy += (Data[0][i])*(Data[1][i]);
        sum_x_square += (Data[0][i])*(Data[0][i]);
    }

    //Creating and then filling two vectors for the q and p values to then fill the 2d vector. this is probably a dumb way to do this
    std::vector<double> coeffs;
    std::vector<double> chi_square_vec;

    //Calculating p and q from eqt.
    p = (sum_xy - sum_x*sum_y)/(sum_x_square - sum_x*sum_x);
    q = (sum_x_square*sum_y - sum_xy*sum_x)/(Data[0].size()*sum_x_square - sum_x*sum_x);

    //Creating a string of the function to write to a file:
    std::string p_string = std::to_string(p);
    std::string q_string = std::to_string(q);
    std::string function = "y = " + p_string + "x" + " + " + q_string;

    //The two coeffecients for the linear fit are added to the vector[0][] position:
    coeffs.push_back(p);
    coeffs.push_back(q);

    best_fits.push_back(coeffs);

    //Chi_squared analysis:
    //Initialising chi_squared variable:
    double chi_squared;

    for(int i=0; i<Data[0].size(); i++){
        auto observed = Data[1][i];
        auto expected = p*Data[0][i] + q;
        auto meas_error = error[1][i];  //Error on y
        chi_squared += ((observed-expected)*(observed-expected))/(meas_error*meas_error);
        //std::cout << chi_squared << std::endl;   //Print statement for debugging
    }

    //The chi squared value is added to it's own vector which is pushed back to the vector[1][] position
    chi_square_vec.push_back(chi_squared);
    best_fits.push_back(chi_square_vec);


    return best_fits;

}


//Fucntion below takes a set of (x,y) data as an input and returns x^y (without using a for loop of pow() fucntion)
std::vector<double> x_pow_y(std::vector<std::vector<double>> data){
    //Intiialising a variable equal to the number of data points:
    int i = data[0].size();

    //Initialising the vector to be returned which will hold the x^y vals:
    std::vector<double> x_pow_y;

    while(i>0){
        //Need to create a index variable that is the difference between the number of data points and the current i value, otherwise operation is ran in reverse order:
        int index = data[0].size() - i;

        double x = data[0][index];
        double x_temp = data[0][index];
        int y = int(data[1][index]);   //This always rounds down: You need to find a way to do rounding correctly
        //std::cout << "x = " << x << ", y = " << y << std::endl;  //Print statement for debugging

        //Looping through each value of y and multiplying x by itself:
        while(y>1){
            x_temp = x_temp*x;
            y--;
        }
        x_pow_y.push_back(x_temp);
        i--;
    }

    return x_pow_y;
}


//Function which takes a user selection input and returns the result:
int print_save_selection(){
    //Initialising selection variables:
    int selection;

    std::cout << "Please enter your selection: " << "\n" << "[1] Print" << "\n" << "[2] Save to file" << std::endl;
    std::cin >> selection;

    if(selection < 0 || selection > 2){
        std::cout << "Invalid input" << std::endl;
        return 0;
    }
    else{
        return selection;
    }
}


//Functions to save data to files:

//For saving function/chi squared
void save(std::vector<std::vector<double>> data, std::string header){
    //Extracting the p, q and chi squared variables
    auto p = data[0][0];
    auto q = data[0][1];
    auto chi_squared = data[1][0];

    std::string function("y = " + std::to_string(p) + "x + " + std::to_string(q));

    //Code below creates a new text file and adds the function y = px + q and chi squared values:
    std::ofstream newFile("function_and_chi_squared.txt");
    //Adding the header:
    newFile << header << "\n";
    newFile << "Function: " << function << "\n" << "Chi squared = " << std::to_string(chi_squared);
    newFile.close();
}

//For saving magnitude/x^y
void save(std::vector<double> data, std::string header, std::string fileName){
    //Creating new file:
    std::ofstream newFile(fileName);
    newFile << header << "\n";

    for(int i=0; i<data.size(); i++){
        newFile << data[i] << "\n";
    }

    newFile.close();
}