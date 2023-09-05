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

#include <vector>
#include <array>
#include <boost/multi_array.hpp>

struct axis
{
    uint32_t nx;
    double   dx;
    double   idx;
    std::vector<double> x;

    void reset()
    {
        nx = 0;
        dx = 0;
        idx = 0;
        x.clear();
    }
};

typedef std::array<axis, 2>           axes2D;
typedef boost::multi_array<double, 2>        data2D;

