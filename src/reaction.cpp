/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "reaction.h"

#include <cmath>
#include <iostream>
#include <plog/Log.h>

double
reaction::rate(double tgas, double cosmic)
{
  int type = int(a_rad);
  double k = 0.0;
  switch (type) {
    case 1:
      k = alpha * cosmic;
      break;
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 14:
      k = alpha * pow(tgas / 300.0, beta) *
          exp(-sigma / tgas); //-sigma was -gamma, this should only be called
                              //for chem reactions
      break;
    default:
      PLOGE << "UNKNOWN TYPE (" << type << ")";
  }

  return k;
}
