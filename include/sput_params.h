/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#pragma once

#include "sputter.h"

#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <fstream>

struct params
{
    std::vector<std::vector<double>> eiCoeff;
    std::vector<std::vector<double>> Eth;
    std::vector<std::vector<double>> asc;
    std::vector<std::vector<double>> SiCoeff;
    std::vector<std::vector<double>> mu;
    std::vector<std::vector<double>> alpha;
    std::vector<double> mi;
    std::vector<double> zi;
    std::vector<double> miGRAMS;
    std::vector<double> miKG;
    std::vector<double> md;
    std::vector<double> mdGRAMS;
    std::vector<double> rhod;
    std::vector<double> zd;
    std::vector<double> u0;
    std::vector<double> K;
    std::vector<double> msp_2rhod;
    std::vector<double> y8_piMi; 
    std::vector<double> three_2Rhod;

    void alloc_vecs(std::size_t numGrns, std::size_t numGas)
    {
        mi.resize(numGas);
        y8_piMi.resize(numGas);
        zi.resize(numGas);
        miGRAMS.resize(numGas);
        miKG.resize(numGas);
        md.resize(numGrns);
        mdGRAMS.resize(numGrns);
        rhod.resize(numGrns);
        zd.resize(numGrns);
        u0.resize(numGrns);
        K.resize(numGrns);
        msp_2rhod.resize(numGrns);
        three_2Rhod.resize(numGrns);

        mu.resize(numGrns);
        eiCoeff.resize(numGrns);
        Eth.resize(numGrns);
        asc.resize(numGrns);
        alpha.resize(numGrns);
        SiCoeff.resize(numGrns);

        std::fill(Eth.begin(), Eth.end(), std::vector<double>(numGas, 0.0));

        for (size_t i = 0; i < numGrns; ++i)
        {    
            mu[i].resize(numGas);
            eiCoeff[i].resize(numGas);
            asc[i].resize(numGas);
            alpha[i].resize(numGas);
            SiCoeff[i].resize(numGas);
        }
    }

};
