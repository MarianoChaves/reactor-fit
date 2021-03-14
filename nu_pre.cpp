#include "nu_pre.h"
#include <stdio.h>
#include <stdlib.h>


Pre::Pre(double(*f)(double arg)){
    this->setFunction(f);
}

Pre::Pre(double(*f)(double arg), std::vector<double> q){
    this->setFunction(f);
    this->setQ(q);
}
Pre::~Pre(){}

void Pre::calc(int N)
{
    F.clear();
    for(int i=0; i<q_size; i++)
    {
        F.push_back(function(q[i]));
    }
}

void Pre::setGeo(double(*geometry)(double(*f)(double arg), double arg))
{
    F.clear();
    for(int i=0; i<q_size; i++)
    {
        F.push_back(geometry(function,q[i]));
    }
}
