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

#include "reaction.h"

//#define BOOST_SPIRIT_DEBUG
#include <map>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <boost/fusion/adapted.hpp>
#include <boost/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <cmath>

// #include "bilinear_interpolator.h"
// #include "bicubic_interpolator.h"

namespace qi  = boost::spirit::qi;
namespace phx = boost::phoenix;

typedef std::vector<reaction> reaction_v;

BOOST_FUSION_ADAPT_STRUCT(
  reaction,
  (spec_v, reacts)(spec_v, prods)(spec_v, ks_list)(double, alpha)(double, beta)(
    double,
    sigma)(double, a_rad)(int, type)(std::string, extra))

template<typename Iterator, typename Skipper>
struct network_parser : qi::grammar<Iterator, reaction_v(), Skipper>
{
  network_parser() : network_parser::base_type(network_rule)
  {
    using qi::char_;
    using qi::double_;
    using qi::eoi;
    using qi::graph;
    using qi::int_;
    using qi::lexeme;
    using qi::lit;
    using qi::omit;
    using qi::space;

    spec_rule     = lexeme[+graph] % '+';
    extra_rule    = +(char_ - '\n');
    reaction_rule = spec_rule >> lit("->") >> spec_rule >> lit("|") >>
                    spec_rule >> double_ >> double_ >> double_ >> double_ >>
                    int_ >> *(extra_rule);

    network_rule = reaction_rule % +char_("\n") >> omit[*space] > eoi;

    BOOST_SPIRIT_DEBUG_NODES(
      (network_rule)(reaction_rule)(spec_rule)(extra_rule));
  }
  qi::rule<Iterator, reaction_v(), Skipper> network_rule;
  qi::rule<Iterator, reaction(), Skipper> reaction_rule;
  qi::rule<Iterator, spec_v(), Skipper> spec_rule;
  qi::rule<Iterator, std::string(), Skipper> extra_rule;
};

// typedef bilinear_interpolator interpolator;
// typedef bicubic_interpolator interpolator;


struct network
{
  std::vector<reaction> reactions;
  std::vector<std::string> species;
  std::vector<std::vector<size_t>> reactants_idx;
  std::vector<std::vector<size_t>> products_idx;
  std::vector<std::vector<size_t>> ks_lists_idx;
  std::vector<size_t> chemical_reactions_idx;
  std::vector<size_t> nucleation_reactions_idx;
  std::vector<size_t> gases_idx;
  std::vector<std::map<size_t, uint32_t>> nucleation_species_count;

  // std::map<int, interpolator> nucl_rate_data;
  std::string network_label;
  size_t n_species = 0;
  size_t n_reactions = 0;
  size_t n_nucleation_reactions;
  size_t n_chemical_reactions;

  // sms added this here, may move
  std::vector<std::string> grn_names;


  void get_species_list();
  void map_species_to_reactions();
  int get_species_index(const std::string& spec) const;
  void read_network(const std::string& chemfile);
  void post_process();
  network();
  virtual ~network();
};
