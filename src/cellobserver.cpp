/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "cellobserver.h"

#include "network.h"
#include "cell.h"

#include <plog/Log.h>
#include <string>
#include <vector>
#include <sstream>
#include <boost/filesystem.hpp>

// intialize the writer class
CellObserver::CellObserver(std::size_t cid,const network* net, configuration* con)
  : cid(cid)
{
  num_nuc  = net->n_nucleation_reactions;
  num_spec = net->n_species;
  numBins = con->bin_number;
  modNum = con->mod_number; 
  m_nstore = con->io_disk_n_steps;
  m_ndump = con->io_screen_n_steps;
  ofname = "output/mod_" + modNum+"/cmNoCool_B"+std::to_string(numBins)+"_m"+modNum+"_"+net->network_label +"_"+std::to_string(cid)+".dat";    
  RSname = "restart/mod_" + modNum+"/restartNC_B"+std::to_string(numBins)+"_m"+modNum+"_"+net->network_label +"_"+std::to_string(cid)+".dat";   
}

// write data to file
void
CellObserver::dump_data()
{
  ofs.open(ofname, std::ofstream::app);
  boost::format fmt("%1$5e ");


  boost::format fmtL("%1$9e ");

  ofs << fmtL % m_state.time << "\n";
  for(double &val: m_state.sizeBins)
    ofs << fmtL % val << " ";
  ofs << "\n";

  ofs.close();
}

// write data to restart file
void
CellObserver::restart_dump ()
{
  oRS.open(RSname);
  boost::format fmt("%1$5e ");
  boost::format fmtL("%1$9e ");

  oRS << fmtL % m_state.time << "\n";
  for(double &val: m_state.init_abund)
    oRS << fmtL % val << " ";
  oRS << "\n";
  for(double &val: m_state.abund_moments_sizebins)
    oRS << fmtL % val << " ";
  oRS << "\n";
  for(double &val: m_state.vd)
    oRS << fmtL % val << " ";
  oRS << "\n";
  for(double &val: m_state.sizeBins)
    oRS << fmtL % val << " ";
  oRS << "\n";
  for(double &val: m_state.erosion_dadt)
    oRS << fmtL % val << " ";
  oRS << "\n";
  for(double &val: m_state.runningTot_size_change)
    oRS << fmtL % val << " ";
  oRS << "\n";

  oRS.close();
}

// call to class, if user specified, write to file or restart file
void
CellObserver::operator()(const cell_state& s)
{
  ++n_called;
  if (n_called % m_nstore == 0) {
    m_state = s;
    dump_data();
    restart_dump();
  }
}

// final writing of data
void CellObserver::finalSave (const cell_state& s)
{   
  dump_data();
  boost::filesystem::remove(RSname);
}