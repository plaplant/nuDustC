#include <fstream>
#include <vector>
#include <unordered_set>
#include <string>
#include <cmath>

#include <boost/filesystem/path.hpp>
#include <boost/spirit/include/qi.hpp>
#include <plog/Log.h>

#include "network.h"


double M_Pi = 3.141592;

network::network() :
    n_reactions ( 0 ),
    n_species ( 0 )
{
}

network::~network()
{

}

// read the network file and update the species list
void
network::read_network ( const std::string &chemfile )
{
    //PLOGI << "Reading network from [" << chemfile << "]";

    std::ifstream ifs ( chemfile.c_str(), std::ifstream::in );
    ifs.unsetf ( std::ios::skipws );
    boost::spirit::istream_iterator f ( ifs ), l;

    //PLOGI << "Network file [" << chemfile << "] opened";

    network_parser<boost::spirit::istream_iterator, qi::blank_type> p;

    //PLOGI << "parser applied";

    bool ok = qi::phrase_parse ( f, l, p, qi::blank, reactions );

    if ( ok )
    {
        //PLOGI << "parser ok";
    }
    else
    {
        PLOGE << chemfile << " parser fail!!";
        exit ( 1 );
    }

    //if ( f != l ) PLOGI << "remaining unparsed: '" << std::string ( f, l ) << "'";

    //PLOGI << "getting species from network...";
    get_species_list();

    n_reactions = reactions.size();
    n_species = species.size();

    //PLOGI << "mapping species index to reactions...";
    map_species_to_reactions();

    //PLOGI << "network loaded from [" << chemfile << "]";
    //PLOGI << "got " << n_species << " species in " << n_reactions << " reactions";

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
    for ( auto i = 0 ; i < n_reactions; ++i )
    {
        if ( reactions[i].type == REACTION_TYPE_NUCLEATE )
        {
            //PLOGI << "reaction " << reactions[i].id << " is dust nucleation";
            
            std :: unordered_set<size_t> react_elems;

            for ( auto &re : reactants_idx[i] ) react_elems.emplace ( re );
           
            std :: map < size_t, uint32_t > tmp_map;

            for ( auto &re : react_elems )
            {
              tmp_map.insert( std::make_pair ( re,  std :: count (  std::begin (reactants_idx[i]),
                                                                    std::end   (reactants_idx[i]),
                                                                    re ) ) );
            }

            nucleation_species_count.emplace_back ( tmp_map );
            nucleation_reactions_idx.emplace_back ( i );
            //PLOGI << "nucleation data loaded";
            /*
            for ( auto &kv : tmp_map)
            {
              PLOGI << "species " << kv.first << " has " << kv.second << " count";
            }
            */

            double ang2cm = 1.0E-8;
            reactions[i].a_rad = reactions[i].a_rad*ang2cm;
            reactions[i].omega0 = pow(reactions[i].a_rad,3) * 4.0 / 3.0 * M_Pi;
            //reactions[i].alpha = reactions[i].alpha * 1E4;
        }
        else
        {
            chemical_reactions_idx.emplace_back ( i );
        }
    }
    n_nucleation_reactions = nucleation_reactions_idx.size();
    n_chemical_reactions = chemical_reactions_idx.size();
    //PLOGI << "Nucleation Reactions: "<<n_nucleation_reactions<<" Chemical Reactions: "<<n_chemical_reactions;
}

/*
 * returns the internal index that corrisponds to the
 * species
 *
 */
int
network::get_species_index ( const std::string &spec ) const
{
    auto it = std::find ( std :: begin(species), std :: end(species), spec );
    auto idx = -1;
    if ( it != species.end() )
        idx = std::distance ( std :: begin(species), it );
    //else
        //PLOGE << "ERROR: cannot find species(" << spec << ")";
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
    reactants_idx.resize ( n_reactions );
    products_idx.resize ( n_reactions );
    ks_lists_idx.resize ( n_reactions );
    for ( auto i = 0; i < n_reactions; ++i )
    {
        for ( const auto &reactant : reactions[i].reacts )
        {
            reactants_idx[i].push_back ( get_species_index ( reactant ) );
            if (std::find(gases_idx.begin(), gases_idx.end(), get_species_index ( reactant )) == gases_idx.end())
                gases_idx.push_back(get_species_index ( reactant ));
        }
        for ( const auto &product : reactions[i].prods )
        {
            products_idx[i].push_back ( get_species_index ( product ) );
            if ( reactions[i].type == REACTION_TYPE_NUCLEATE )
                grn_names.push_back(product);
        }

        for ( const auto &k_spec : reactions[i].ks_list )
            ks_lists_idx[i].push_back ( get_species_index ( k_spec ) );
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

    for ( const auto &r : reactions )
    {
        spec_set.insert ( r.reacts.begin(), r.reacts.end() );
        spec_set.insert ( r.prods.begin(), r.prods.end() );
    }
    //spec_v oxy{"CO","SiO"};
    //spec_set.insert(oxy.begin(),oxy.end());
    species.insert ( species.begin(), spec_set.begin(), spec_set.end() );

    auto s_id = 0;

}

