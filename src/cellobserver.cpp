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
  ofname = "output/B"+std::to_string(numBins)+"_"+net->network_label +"_"+std::to_string(cid)+".dat";    
  RSname = "restart/restart_B"+std::to_string(numBins)+"_"+net->network_label +"_"+std::to_string(cid)+".dat"; 

  for(auto gn_id =0; gn_id < num_nuc; gn_id++)
  {
      grnNames += net->reactions[gn_id].prods[0] +" ";
  }  
}

void CellObserver::init_dump(const cell_state& s)
{
    m_state = s;
    ofs.open(ofname);
    ofs << grnNames << "\n";

    boost::format fmt ( "%1$14e " );
    boost::format fmtL("%1$9e ");

    for(double &val: m_state.grn_sizes)
    {
      ofs << fmt % val << " ";
    }
    ofs << "\n";
    for(double &val: m_state.abund_moments_sizebins)
    {
      ofs << fmtL % val << " ";
    }
    ofs << "\n";
    ofs.close();
}

// write data to file
void
CellObserver::dump_data(const cell_state& s)
{
  m_state = s;
  ofs.open(ofname, std::ofstream::app);
  boost::format fmt("%1$5e ");
  boost::format fmtL("%1$9e ");

  ofs << fmtL % m_state.time << "\n";
  for(double &val: m_state.abund_moments_sizebins)
    ofs << fmtL % val << " ";
  ofs << "\n";

  ofs.close();
}

// write data to restart file
void
CellObserver::restart_dump(const cell_state& s)
{
    oRS.open(RSname);
    boost::format fmt("%1$5e ");
    boost::format fmtL("%1$9e ");
    oRS << (m_state.time) << "\n";

    for(int idx = 0; idx < m_state.vd.size(); idx ++)
    {
        oRS << fmt % (m_state.vd[idx]) << " ";
    }
    oRS << "\n";

    for(int idx = 0; idx < m_state.runningTot_size_change.size(); idx ++)
    {
        oRS << fmt % (m_state.runningTot_size_change[idx]) << " ";
    }
    oRS << "\n";
    // solution vector
    for(int idx = 0; idx < m_state.abund_moments_sizebins.size(); idx ++)
    {
        oRS << fmt % (m_state.abund_moments_sizebins[idx]) << " ";
    }
    oRS << "\n";

    oRS.close();
}

// call to class, if user specified, write to file or restart file
void
CellObserver::operator()(const cell_state& s)
{
  ++n_called;
  if (n_called % m_ndump == 0) {
    m_state = s;
    dump_data(s);
    restart_dump(s);
  }
}

// final writing of data
void CellObserver::finalSave (const cell_state& s)
{   
  m_state = s;
  dump_data(s);
  boost::filesystem::remove(RSname);
}
