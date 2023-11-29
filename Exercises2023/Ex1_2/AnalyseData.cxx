//Author: Rory McFeely
//Date: 15/11/23
//This code potentially does something

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include "CustomFunctions.h"



int main(){

    std::string line;
    std::string fileName = "input2D_float.txt";
    std::string errorFileName = "error2D_float.txt";
    
    auto data = extract(fileName);
    auto error = extract(errorFileName);
    auto fit_vals = fit_function(data, error);

    for(int i=1; i>0;){
    //Initialising variable for user selection
    int selection;

    //Taking user input to determine which function to perform:
    std::cout << "Welcome to the salty spitoon what'll it be?" << "\n" << "[1] Print raw data" << "\n" << "[2] Calculate magnitudes" << "\n" << "[3] Calculate best fit coeffiecients and chi squared" << "\n" << "[4] Calculate x^y" << std::endl;
    std::cin >> selection;


    //This is probably better done with a switch statement, but didn't have time to implement
    if (selection < 0 || selection > 4){
        std::cout << "Invalid input" << std::endl;
    }
    else if(selection == 1){
        print(data);
    }
    else if(selection == 2){
        auto magnitudes = magnitude(data);
        int print_save = print_save_selection();
        if(print_save == 1){
            print(magnitudes, std::string("Magntidues:"));
        }
        else if(print_save == 2){
            save(magnitudes, std::string("Calculated magnitudes"), std::string("magnitudes.txt"));
        }

    }
    else if(selection == 3){    
        auto fit_results = fit_function(data, error);
        int print_save = print_save_selection();
        if(print_save == 1){
            print(fit_results, std::string("Results of fit function:"));
        }
        else if(print_save == 2){
            save(fit_results, std::string("Results of fit"));
        }          
    }
    else if(selection == 4){
        auto x_topow_y = x_pow_y(data);
        int print_save = print_save_selection();
        if(print_save == 1){
            print(x_topow_y, std::string("x^y:"));
        }
        else if(print_save == 2){
            save(x_topow_y, std::string("x^y"), std::string("x_to_power_y.txt"));
        }
    }

    //The lines below take a user input to determine whether to break from the loop or perform another command
    std::string user_selection;
    std::cout << "Would you like to perform another action? [y/n]" << std::endl;
    std::cin >> user_selection;
    if(user_selection == "n"){
        break;
    }

    }

    return 0;
}

