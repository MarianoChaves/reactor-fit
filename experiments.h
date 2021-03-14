#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include<math.h>
#include<stdlib.h>
#include<time.h>

/*
    This file contain various namespaces 
    that can be useful for neutrino physics.

    --geo: gives geometry of several experiments
*/

namespace geo{
    double mean3D(double(*f)(double arg), double L, double r, double q){
        double mean = 0;
        int sample = 200;
        std::srand(time(NULL));
        for(int i=0; i<sample; i++){

            double rai = (std::rand() % 31)*r/30;
            double thi = (std::rand() % 31)*2*M_PI/30;
            double phi = (std::rand() % 31)*M_PI/30;
            double xl = sqrt(L*L+rai*rai-2*L*rai*sin(phi)*sin(thi));
            mean = mean+f(q*xl);
        }
        mean = mean/sample;
        return mean;
    }

    double pontual_reactor_set_DB(double *Ld, double *td, int Nr, 
    double(*f)(double arg), double q){

        double norm = 0;
        double F = 0;

        double P[3][6] = {{2.082, 2.874, 2.516, 2.554, 2.825, 1.976},
            {2.082, 2.874, 2.516, 2.554, 2.825, 1.976},
            {2.51, 2.45, 2.57, 2.52, 2.52, 2.55}};

        for(int r = 0; r<Nr; r++){
            for(int s=0; s<3; s++){
                norm += td[s]*P[s][r]/(Ld[r]*Ld[r]);
                F += td[s]*P[s][r]*mean3D(f,Ld[r]/1000,3.0/1000,q)/(Ld[r]*Ld[r]); 
                // by 1000 division is converting  meters to km
            };
        }
        return F/norm;
    }
    
    double pontual_reactor_set_DC(double *Ld, double *Pr, int Nr, 
    double(*f)(double arg), double q){

        double norm = 0;
        double F = 0;


        for(int r = 0; r<Nr; r++){
            norm += Pr[r]/(Ld[r]*Ld[r]);
            F += Pr[r]*mean3D(f,Ld[r]/1000,3.0/1000,q)/(Ld[r]*Ld[r]); 
            // by 1000 division is converting  meters to km
        };

        return F/norm;
    }
    
    double pontual_reactor_set_RN(double *Ld, double *Pr, int Nr, 
    double(*f)(double arg), double q){

        double norm = 0;
        double F = 0;


        for(int r = 0; r<Nr; r++){
            norm += Pr[r]/(Ld[r]*Ld[r]);
            F += Pr[r]*mean3D(f,Ld[r]/1000,3.0/1000,q)/(Ld[r]*Ld[r]); 
            // by 1000 division is converting  meters to km
        };

        return F/norm;
    }

    double pontual_reactor_set_KL(double *Ld, double *Pr, int Nr, 
    double(*f)(double arg), double q){

        double norm = 0;
        double F = 0;
        double Lmean = 0;


        for(int r = 0; r<Nr; r++){
            norm += Pr[r]/(Ld[r]*Ld[r]);
            Lmean += (Pr[r]/(Ld[r]*Ld[r]))*Ld[r];
            F += Pr[r]*mean3D(f,Ld[r]/1000,3.0/1000,q)/(Ld[r]*Ld[r]); 
            // by 1000 division is converting  meters to km
        };
        return F/norm;
    }

/* DAYA BAY */
    /*H1*/
    double DB_AD1_BL[6] = {362.38,371.76,903.47,817.16,1353.62,1265.32};
    double DB_AD2_BL[6] = {357,368,903,816,1354,1265};
    /*H2*/
    double DB_AD3_BL[6] = {1332,1358,467,489,557,499};
    double DB_AD8_BL[6] = {1337.43,1362.88,472.97,495.35,558.71,501.07};
    /*H3*/
    double DB_AD4_BL[6] = {1919,1894,1533,1533,1551,1524};
    double DB_AD5_BL[6] = {1917,1891,1533,1533,1554,1528};
    double DB_AD6_BL[6] = {1925,1899,1538,1539,1556,1530};
    double DB_AD7_BL[6] = {1923,1897,1540,1540,1559,1533};
    
    double t1[3]={217,0,1524};
    double t2[3]={217,217,1524};
    double t3[3]={217,217,1524};
    double t4[3]={217,217,1524};
    double t5[3]={217,217,1524};
    double t6[3]={217,217,1524};
    double t7[3]={0,217,1524};
    double t8[3]={0,217,1524};

    double DB_AD1(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD1_BL, t1, Nr, f, q);
    };
    double DB_AD2(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD2_BL, t2, Nr, f, q);
    };
    double DB_AD3(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD3_BL, t3, Nr, f, q);
    };
    double DB_AD4(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD4_BL, t4, Nr, f, q);
    };
    double DB_AD5(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD5_BL, t5, Nr, f, q);
    };
    double DB_AD6(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD6_BL, t6, Nr, f, q);
    };
    double DB_AD7(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD7_BL, t7, Nr, f, q);
    };
    double DB_AD8(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_DB(DB_AD8_BL, t8, Nr, f, q);
    };
    
    /* Double Chooz */
    double DC_ND_BL[2] = {355, 469};
    double DC_FD_BL[2] = {998, 1115};

    double P_ND[2] = {8.5, 8.5};
    double P_FD[2] = {8.5, 8.5};

    double DC_ND(double(*f)(double arg), double q){
        int Nr = 2;
        return pontual_reactor_set_DC(DC_ND_BL, P_ND, Nr, f, q);
    };
    double DC_FD(double(*f)(double arg), double q){
        int Nr = 2;
        return pontual_reactor_set_DC(DC_FD_BL, P_FD, Nr, f, q);
    };

/* RENO */
    double RENO_ND_BL[6] = {660.064, 444.727, 301.559, 339.262, 519.969, 746.155};
    double RENO_FD_BL[6] = {1563.771, 1460.826, 1397.813, 1380.062, 1409.389, 1483.001};

    double P_ND_RN[6] = {2.381, 2.084, 2.224, 2.060, 2.315, 2.248};
    double P_FD_RN[6] = {2.363, 2.082, 2.144, 2.111, 2.387, 2.319};

    double RENO_ND(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_RN(RENO_ND_BL, P_ND_RN, Nr, f, q);
    };
    double RENO_FD(double(*f)(double arg), double q){
        int Nr = 6;
        return pontual_reactor_set_RN(RENO_FD_BL, P_FD_RN, Nr, f, q);
    };

/* KAMLAND */
	double KAMLAND_BL[20]={160000.0, 179500.0, 190600.0, 214000.0, 138600.0, 80600.0, 145400.0, 344000.0, 344000.0, 294600.0, 414000.0, 430200.0, 561200.0, 755400.0, 824100.0, 783500.0, 750000.0, 690000.0, 940000.0, 700000.0};
    double P_KL[20]={24.6,13.7,10.2,10.6,4.5,1.6,4.9,14.2,13.2,3.3,3.8,4.8,6.0,6.7,3.3,5.3,11.2,8.1,16.8,8.9};

	//double KAMLAND_BL[20]={160000,179000,191000,138000,214000,146000,88000,349000,345000,295000,431000,401000,561000,755000,830000,783000,712000,986000,735000,709000};
	//double P_KL[20]={24.3,13.7,10.2,4.5,10.6,4.9,1.6,14.2,13.2,3.3,6.5,3.8,6.0,10.1,5.3,3.3,11.5,17.4,9.2,8.2};
    //double KAMLAND_BL[20]={180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000,180000};


    double KAMLAND(double(*f)(double arg), double q){
        int Nr = 20;
        return pontual_reactor_set_KL(KAMLAND_BL, P_KL, Nr, f, q);
    };

/* */

/* */
}


#endif
