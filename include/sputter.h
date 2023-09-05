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

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <typeinfo> 

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <plog/Log.h>
#include <boost/foreach.hpp>

namespace sput
{

struct GRN_t
{
  std::string name;
  double u0, md, zd, K, mu, rhod;
  explicit GRN_t(std::string _name, double _u0, double _md, double _zd, 
               double _K, double _mu, double _rhod )
      : name(std::move( _name )), u0(_u0), md(_md), zd(_zd), K(_K), mu(_mu), rhod(_rhod)
  { }
};

struct grain_comp_t
{
  std::string name;
  std::vector<int> reactAMT;
  std::vector<std::string> react;
  explicit grain_comp_t(std::string _name, std::vector<int> _reactAMT, std::vector<std::string> _react)
      : name(std::move( _name )), reactAMT(_reactAMT), react(std::move( _react ))
  { }
};

struct ions_t
{
  std::string name;
  double zi;
  double mi;
  explicit ions_t(std::string _name, double _zi, double _mi)
      : name(std::move( _name )), zi(_zi), mi(_mi)
  { }
};

struct sputter_list_t
{
  std::map<std::string, GRN_t> grains;
  std::map<std::string, grain_comp_t> grn_comp;
  std::map<std::string, ions_t> ions;
  sputter_list_t(const std::string& sputter_file )
  {
    namespace pt = boost::property_tree;
    pt::ptree root;
    pt::read_json( sputter_file, root );
    for ( pt::ptree::value_type& gr : root.get_child( "data" ) )
    {
      auto rec = gr.second;
      grains.try_emplace(
          rec.get<std::string>( "name" ),
          rec.get<std::string>( "name" ),
          rec.get<double>("u0", -1.f), rec.get<double>("md", -1.f),
          rec.get<double>("zd", -1.f), rec.get<double>("K", -1.f),
          rec.get<double>("mu", -1.f), rec.get<double>("rhod", -1.f));
    }
    for ( pt::ptree::value_type& gr : root.get_child( "grainsCOMP" ) )
    {
      std::vector<int> reactAMT ={};
      std::vector<std::string> react ={};
      BOOST_FOREACH(boost::property_tree::ptree::value_type &amt, gr.second.get_child("reacAMT."))
      {
        reactAMT.push_back(std::stoi(amt.second.data()));
      }
      BOOST_FOREACH(boost::property_tree::ptree::value_type &coeff, gr.second.get_child("react."))
      {
        std::string temp = coeff.second.data();
        react.push_back(temp);
      }

      // todo: add react elem string vector like what was done with reacAMT
      
      auto rec = gr.second;
      grn_comp.try_emplace(
          rec.get<std::string>( "name" ),
          rec.get<std::string>( "name" ),
          reactAMT,react);
    }
    for ( pt::ptree::value_type& gr : root.get_child( "ions" ) )
    {
      auto rec = gr.second;
      ions.try_emplace(
          rec.get<std::string>( "name" ),
          rec.get<std::string>( "name" ),
          rec.get<double>("zi"),
          rec.get<double>("mi"));
    }
  }

  int
  get_grn_index(const std::string& en )
  {
    auto eitr = grains.find( en );
    return (eitr == grains.end() ? -1 : std::distance( grains.begin(), eitr ) );
  }
  int
  get_grnCOMP_index(const std::string& en )
  {
    auto eitr = grn_comp.find( en );
    return (eitr == grn_comp.end() ? -1 : std::distance( grn_comp.begin(), eitr ) );
  }
  int
  get_ions_index(const std::string& en )
  {
    auto eitr = ions.find( en );
    return (eitr == ions.end() ? -1 : std::distance( ions.begin(), eitr ) );
  }
};

}
