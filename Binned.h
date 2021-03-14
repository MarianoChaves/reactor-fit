#ifndef BINNED_H
#define BINNED_H

#include<string>
#include<stdio.h>
#include<vector>
#include<gsl/gsl_spline.h>
#include<math.h>

class Binned{
        std::vector<double> true_value;
        std::vector<double> true_energy;
    public:
        std::vector<double> value;
        std::vector<double> energy;

        void setTrues(){true_value=value; true_energy=energy;};
        void interpolate(std::vector<double> energy);

        Binned(double ini, double fin, int N, std::string spacing);
        Binned(double ini, double fin, int N);
        Binned(){};
        ~Binned(){};
        
};

class input: public Binned{
    private:
        std::vector<double> true_value;
        std::vector<double> true_energy;
        
        int x_collum = 1;
        int y_collum = 2;

        void read_input(char *file_name);
    public:
        char header[12] = "x,y\n";

        input(char *file_name, int x, int y);
        ~input(){};
};




#endif