/*© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.*/

#include "network.h"

#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/qi.hpp>
#include <cmath>
#include <fstream>
#include <plog/Log.h>
#include <string>
#include <unordered_set>
#include <vector>


double M_Pi = 3.141592;

network::network() :
    n_reactions ( 0 ),
    n_species ( 0 )
{
}

network::~network() {}

void
network::read_network(const std::string& chemfile)
{
  std::ifstream ifs(chemfile.c_str(), std::ifstream::in);
  ifs.unsetf(std::ios::skipws);
  boost::spirit::istream_iterator f(ifs), l;

  network_parser<boost::spirit::istream_iterator, qi::blank_type> p;

  bool ok = qi::phrase_parse(f, l, p, qi::blank, reactions);
  if (ok) {
    //PLOGI << "network parser ok";
  } else {
    PLOGE << "network parser fail!!";
    exit(1);
  }

  //PLOGI << "getting species from network...";
  get_species_list();

  n_reactions = reactions.size();
  n_species   = species.size();
  map_species_to_reactions();
  network_label = boost::filesystem::path(chemfile).stem().string();
}

/*
 * do any updates on read network
 * used to:
 * - read nucleation interpolators for dust grains
 *
 */
void
network::post_process()
{
  for (size_t i = 0; i < n_reactions; ++i) {
    if (reactions[i].type == REACTION_TYPE_NUCLEATE) {
      // PLOGI << "reaction " << reactions[i].id << " is dust nucleation";

      std ::unordered_set<size_t> react_elems;

      for (auto& re: reactants_idx[i])
        react_elems.emplace(re);

      std ::map<size_t, uint32_t> tmp_map;

      for (auto& re: react_elems) {
        tmp_map.insert(std::make_pair(re,
                                      std ::count(std::begin(reactants_idx[i]),
                                                  std::end(reactants_idx[i]),
                                                  re)));
      }

      nucleation_species_count.emplace_back(tmp_map);
      nucleation_reactions_idx.emplace_back(i);
      // PLOGI << "nucleation data loaded";
      /*
      for ( auto &kv : tmp_map)
      {
        PLOGI << "species " << kv.first << " has " << kv.second << " count";
      }
      */

      double ang2cm       = 1.0E-8;
      reactions[i].a_rad  = reactions[i].a_rad * ang2cm;
      reactions[i].omega0 = pow(reactions[i].a_rad, 3) * 4.0 / 3.0 * M_Pi;
      // reactions[i].alpha = reactions[i].alpha * 1E4;
    } else {
      chemical_reactions_idx.emplace_back(i);
    }
  }
  n_nucleation_reactions = nucleation_reactions_idx.size();
  n_chemical_reactions   = chemical_reactions_idx.size();
  PLOGI << "Nucleation Reactions: " << n_nucleation_reactions
        << " Chemical Reactions: " << n_chemical_reactions;
}

/*
 * returns the internal index that corrisponds to the
 * species
 *
 */
int
network::get_species_index(const std::string& spec) const
{
  auto it  = std::find(std ::begin(species), std ::end(species), spec);
  auto idx = -1;
  if (it != species.end())
    idx = std::distance(std ::begin(species), it);
  return idx;
}

/*
 * maps reactions/products species to their respective
 * indices. this is functionally a std::map, but std::map
 * is slow on lookups, so I decided to code this explicitly.
 */
void
network::map_species_to_reactions()
{
  reactants_idx.resize(n_reactions);
  products_idx.resize(n_reactions);
  ks_lists_idx.resize(n_reactions);
  for (size_t i = 0; i < n_reactions; ++i) {
    for (const auto& reactant: reactions[i].reacts) {
      reactants_idx[i].push_back(get_species_index(reactant));
      if (std::find(gases_idx.begin(),
                    gases_idx.end(),
                    get_species_index(reactant)) == gases_idx.end())
        gases_idx.push_back(get_species_index(reactant));
    }
    for (const auto& product: reactions[i].prods) {
      products_idx[i].push_back(get_species_index(product));
      grn_names.push_back(product);
    }

    for (const auto& k_spec: reactions[i].ks_list)
      ks_lists_idx[i].push_back(get_species_index(k_spec));
  }
}

/*
 * constructs the internal species list from the network.
 * this is a slight hack of std::unordered_set, which
 * only allows unique inserts.
 */
void
network::get_species_list()
{
  std::unordered_set<std::string> spec_set;

  species.clear();

  for (const auto& r: reactions) {
    spec_set.insert(r.reacts.begin(), r.reacts.end());
    spec_set.insert(r.prods.begin(), r.prods.end());
  }
  // spec_v oxy{"CO","SiO"};
  // spec_set.insert(oxy.begin(),oxy.end());
  species.insert(species.begin(), spec_set.begin(), spec_set.end());

  //for (const auto& s: species) {
  //  PLOGI << "added..." << s << " to species list";
  //}
}
