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

#include "H5Cpp.h"
#include "cell.h"

#include <vector>

class CellObserver
{
  std::size_t cid;
  

  uint32_t num_nuc;
  uint32_t num_spec;
  uint32_t numBins;
  std::string modNum;
  std::string   ofname;
  std::string   RSname;

private:

  cell_state m_state;
  double m_time;
  double m_elapsed;

  uint32_t m_nstore, m_ndump;
  uint32_t n_called;

  std::ofstream ofs;
  std::ofstream oRS;
  std::string hfname;

public:
  CellObserver(std::size_t cid, const network* net, configuration* con);
  virtual ~CellObserver() {}

  void operator()(const cell_state& s);
  void dump_data();
  void restart_dump();
  void finalSave(const cell_state& s);
};
