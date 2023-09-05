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

#include <boost/math/interpolators/makima.hpp>


// todo: I could probably get rid of this and declare the boost splines in
// cell.h and cell.cpp
namespace boost::math::interpolators {

template<class RandomAccessContainer>
class makima
{
public:
  using Real = RandomAccessContainer::value_type;

  makima(
    RandomAccessContainer&& abscissas,
    RandomAccessContainer&& ordinates,
    Real left_endpoint_derivative  = std::numeric_limits<Real>::quiet_NaN(),
    Real right_endpoint_derivative = std::numeric_limits<Real>::quiet_NaN());

  Real operator()(Real x) const;

  Real prime(Real x) const;

  void push_back(Real x, Real y);

  friend std::ostream& operator<<(std::ostream& os, const makima& m);
};

} // namespace boost::math::interpolators