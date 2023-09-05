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

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <map>
#include <plog/Log.h>
#include <vector>

namespace xkin {

struct element_t
{
  std::string symbol;
  std::string name;
  std::string category;
  int atomic_number;
  int period;
  double density;
  double mass;
  std::vector<unsigned int> elec_shells;
  explicit element_t(std::string _symbol,
                     std::string _name,
                     std::string _category,
                     int _atomic_number,
                     int _period,
                     double _density,
                     double _mass,
                     std::vector<unsigned int> _elec_shells)
    : symbol(std::move(_symbol)), name(std::move(_name)),
      category(std::move(_category)), atomic_number(_atomic_number),
      period(_period), density(_density), mass(_mass),
      elec_shells(std::move(_elec_shells))
  {
  }
};

struct element_list_t
{
  std::map<std::string, element_t> elements;

  element_list_t(const std::string& elem_file)
  {
    namespace pt = boost::property_tree;

    pt::ptree root;

    pt::read_json(elem_file, root);

    for (pt::ptree::value_type& element: root.get_child("elements")) {
      auto rec = element.second;

      std::vector<unsigned int> shlv;
      for (pt::ptree::value_type& s: rec.get_child("shells"))
        shlv.push_back(std::stoi(s.second.data()));

      elements.try_emplace(
        rec.get<std::string>("symbol"),
        rec.get<std::string>("symbol"), // NB: key, {element stuff}
        rec.get<std::string>("name"),
        rec.get<std::string>("category"),
        rec.get<int>("number"),
        rec.get<int>("period"),
        rec.get<double>("density", -1.f),
        rec.get<double>("atomic_mass", -1.f),
        shlv);
    }
  }

  int get_element_index(const std::string& en)
  {
    auto eitr = elements.find(en);

    return (eitr == elements.end() ? -1
                                   : std::distance(elements.begin(), eitr));
  }
};
} // namespace xkin
