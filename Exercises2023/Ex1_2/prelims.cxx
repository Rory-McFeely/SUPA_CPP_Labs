//Author: Rory McFeely
//Date: 15/11/23
//Preliminary excercise

#include <iostream>
#include <math.h>

double magnitude(double, double);

int main(){
    std::cout << "Hello World!" << std::endl;

    double x = 2.3, y = 4.5;               //Defining x and y variables
    double vectorMag = sqrt((x*x)+(y*y));  //Calculating the magnitude of resultant vector from x and y
    std::cout << "The magnitude of vectors x and y are: " << vectorMag << std::endl;

    auto vectorMag2 = magnitude(x, y);
    std::cout << "The magnitude of vectors x and y using the function are: " << vectorMag2 << std::endl;


    return 0;
}


double magnitude(double x, double y){
    return sqrt((x*x)+(y*y));
}

