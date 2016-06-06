//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_HC_UFO_H
#define HelAmps_HC_UFO_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_HC_UFO 
{
double Sgn(double e, double f); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double>
    fi[18]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void VVS2_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);
void VVS2_6_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex);

void FFS2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);

void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void VVSS1_3(complex<double> V1[], complex<double> V2[], complex<double> S4[],
    complex<double> COUP, double M3, double W3, complex<double> S3[]);
void VVSS1_3_3(complex<double> V1[], complex<double> V2[], complex<double>
    S4[], complex<double> COUP1, complex<double> COUP2, double M3, double W3,
    complex<double> S3[]);

void FFS1_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);
void FFS1_2_1(complex<double> F2[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> F1[]);

void VVS6_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);

void VVVS1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    S4[], complex<double> COUP, double M1, double W1, complex<double> V1[]);
void VVVS1_2P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    S4[], complex<double> COUP1, complex<double> COUP2, double M1, double W1,
    complex<double> V1[]);

void VVSS3_3(complex<double> V1[], complex<double> V2[], complex<double> S4[],
    complex<double> COUP, double M3, double W3, complex<double> S3[]);

void VVV8P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void FFS2_1(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);

void VVVS2P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    S4[], complex<double> COUP, double M1, double W1, complex<double> V1[]);

void VVS6P0_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVS2P0_1(complex<double> V2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);
void VVS2_6P0_1(complex<double> V2[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, double M1, double W1, complex<double> V1[]);

void VVV8_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);

void FFS1_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);
void FFS1_2_2(complex<double> F1[], complex<double> S3[], complex<double>
    COUP1, complex<double> COUP2, double M2, double W2, complex<double> F2[]);

void FFS2_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);

void FFS1_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex);
void FFS1_2_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP1, complex<double> COUP2, complex<double> & vertex);

}  // end namespace MG5_HC_UFO

#endif  // HelAmps_HC_UFO_H
