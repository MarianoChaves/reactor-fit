#ifndef NU_PRE_H
#define NU_PRE_H
#include<vector>
#include<math.h>
/*
    This file of neutrin utilities was modified
    in order to be implemented in GLOBESfit
*/
/*
        This class can be used to calculate
        pre-compilated functions.
        It can depends on the q variables, the 
        geometry of the problem or energy.

        You can calculate f(x) by using set()
        or g(f(x)) by using setGeo(f)
*/

class Pre
{
    private:
        double(*function)(double arg);
        std::vector<double> q;
        int q_size;

    public:
        std::vector<double> F;

        void setFunction(double(*f)(double arg)){function = f;};
        void setQ(std::vector<double> q){this->q=q; q_size = q.size();};

        void calc(int N);

        void setGeo(double(*geometry)(double(*f)(double arg), double arg));
        Pre(double(*f)(double arg), std::vector<double> q);
        Pre(double(*f)(double arg));
        Pre();
        ~Pre();

};

#endif