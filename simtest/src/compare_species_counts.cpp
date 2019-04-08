/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech and Pierre Barbera and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

class ScrappEdgeData : public CommonEdgeData
{
    // -------------------------------------------------------------------
    //     Constructor and Rule of Five
    // -------------------------------------------------------------------

public:

    virtual ~ScrappEdgeData() = default;

    // Move ctor and assignment.
    ScrappEdgeData( ScrappEdgeData&& )             = delete;
    ScrappEdgeData& operator= ( ScrappEdgeData&& ) = delete;

protected:

    ScrappEdgeData() = default;

    // Copy ctor and assignment.
    ScrappEdgeData( ScrappEdgeData const& )             = default;
    ScrappEdgeData& operator= ( ScrappEdgeData const& ) = default;

public:

    static std::unique_ptr< ScrappEdgeData > create()
    {
        return std::unique_ptr< ScrappEdgeData >( new ScrappEdgeData() );
    };

    virtual std::unique_ptr< BaseEdgeData > recreate() const override
    {
        return std::unique_ptr< ScrappEdgeData >( new ScrappEdgeData() );
    }

    virtual std::unique_ptr< BaseEdgeData > clone() const override
    {
        return std::unique_ptr< ScrappEdgeData >( new ScrappEdgeData( *this ));
    }

    // -----------------------------------------------------
    //     Data Members
    // -----------------------------------------------------

    double species_count = 0.0;

};

Tree convert_key_attribute_tree_to_scrapp_tree( AttributeTree const& source )
{
    return convert(
        source,
        [] ( tree::BaseNodeData const& node_data ) {
            auto new_node_data = tree::CommonNodeData::create();
            auto const& orig_node = dynamic_cast< AttributeTreeNodeData const& >( node_data );
            new_node_data->name = orig_node.name;

            return new_node_data;
        },
        [] ( tree::BaseEdgeData const& edge_data ) {
            auto new_edge_data = ScrappEdgeData::create();
            auto const& orig_edge = dynamic_cast< AttributeTreeEdgeData const& >( edge_data );
            new_edge_data->branch_length = orig_edge.branch_length;
            assert( orig_edge.attributes.count( "species_count" ) > 0 );
            new_edge_data->species_count = stod( orig_edge.attributes.at( "species_count" ) );

            return new_edge_data;
        }
    );
}

Tree convert_key_attribute_tree_to_scrapp_mass_tree( AttributeTree const& source )
{
    return convert(
        source,
        [] ( tree::BaseNodeData const& node_data ) {
            auto new_node_data = tree::MassTreeNodeData::create();
            auto const& orig_node = dynamic_cast< AttributeTreeNodeData const& >( node_data );
            new_node_data->name = orig_node.name;

            return new_node_data;
        },
        [] ( tree::BaseEdgeData const& edge_data ) {
            auto new_edge_data = tree::MassTreeEdgeData::create();
            auto const& orig_edge = dynamic_cast< AttributeTreeEdgeData const& >( edge_data );
            new_edge_data->branch_length = orig_edge.branch_length;
            auto bl_half = orig_edge.branch_length / 2.0;
            auto weight = stod( orig_edge.attributes.at( "species_count" ) );
            assert( orig_edge.attributes.count( "species_count" ) > 0 );
            new_edge_data->masses[ bl_half ] = weight;

            return new_edge_data;
        }
    );
}

double calc_per_edge_error( Tree const& lhs, Tree const& rhs )
{
    // need to have identical topology if we're to traverse simultaneously
    if ( not identical_topology( lhs, rhs ) ) {
        throw std::invalid_argument{"Trees have incompatible topology!"};
    }

    std::vector<double> per_edge_error;

    // simultaneous iteration
    auto lhs_it  = postorder( lhs ).begin();
    auto rhs_it  = postorder( rhs ).begin();
    auto lhs_end = postorder( lhs ).end();
    auto rhs_end = postorder( rhs ).end();
    for( ; lhs_it != lhs_end && rhs_it != rhs_end ; ++lhs_it, ++rhs_it ) {

        // If we are at the last iteration, we reached the root. Thus, we have moved all masses
        // and don't need to proceed. If we did, we would count an edge of the root again
        // (because the iterator traverses nodes, not edges, so the root node itself is traversed,
        // although it has no proper edge that we would need to process).
        if( lhs_it.is_last_iteration() || rhs_it.is_last_iteration() ) {
            // If one iterator is at the end, but not the other, something is wrong.
            if( ! lhs_it.is_last_iteration() || ! rhs_it.is_last_iteration() ) {
                throw std::invalid_argument( "Incompatible MassTrees." );
            }
            continue;
        }

        // ================================
        // Actually interesting calculation
        // ================================
        auto lhs_sc = lhs_it.edge().data<ScrappEdgeData>().species_count;
        auto rhs_sc = rhs_it.edge().data<ScrappEdgeData>().species_count;

        per_edge_error.emplace_back( std::abs(lhs_sc - rhs_sc) );
    }

    // Now we need to be done with both trees, otherwise we have a problem.
    if( lhs_it != lhs_end || rhs_it != rhs_end ) {
        throw std::invalid_argument( "Incompatible MassTrees." );
    }

    return std::accumulate( std::begin(per_edge_error), std::end(per_edge_error), 0.0);
}

// =================================================================================================
//     main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    // utils::Logging::log_to_stdout();
    // utils::Logging::details.date = true;
    // utils::Logging::details.time = true;
    // utils::Options::get().number_of_threads( 40 );
    // LOG_BOLD << utils::Options::get().info();
    // LOG_BOLD;

    // Check if the command line contains the right number of arguments.
    if( argc != 3 ) {
      LOG_INFO << "Usage: " << argv[0] << " <scrapp_result (newick)> <true_annotation (newick)>";
      return 1;
    }

    std::string scrapp_tree_file( argv[1] );
    std::string true_tree_file( argv[2] );

    // set up the reader to parse NHX attributes, specifically the ones we want
    auto reader = KeyedAttributeTreeNewickReader();
    reader.set_nhx_delimiters();
    reader.add_attribute( "species_count", KeyedAttributeTreeNewickReaderPlugin::Target::kEdge, "species_count", "0.0" );

    std::vector<MassTree> mass_trees;

    auto attr_tree = reader.read( from_file( scrapp_tree_file ));
    // get two conversions of this tree: once with the values as special species_count
    auto scrapp_tree = convert_key_attribute_tree_to_scrapp_tree( attr_tree );
    // and once already converted into masses on a mass tree, so we can easily get the KR-Distance
    mass_trees.push_back( convert_key_attribute_tree_to_scrapp_mass_tree( attr_tree ) );

    // repeat for the true tree
    attr_tree = reader.read( from_file( true_tree_file ));
    auto true_tree = convert_key_attribute_tree_to_scrapp_tree( attr_tree );
    mass_trees.push_back( convert_key_attribute_tree_to_scrapp_mass_tree( attr_tree ) );

    if ( not identical_topology( mass_trees[0], mass_trees[1] ) ) {
        throw std::invalid_argument{"Trees have incompatible topology!"};
    }

    mass_trees_make_average_branch_lengths( mass_trees );

    // normalize branch lengths
    // scale_all_branch_lengths( mass_trees[0], 1 / length( mass_trees[0] ) );
    // mass_tree_center_masses_on_branches( mass_trees[0] );

    // scale_all_branch_lengths( mass_trees[1], 1 / length( mass_trees[1] ) );
    // mass_tree_center_masses_on_branches( mass_trees[1] );

    std::cout << "pure KRD: " << earth_movers_distance( mass_trees[0], mass_trees[1] ) << std::endl;

    // mass normalization
    mass_tree_normalize_masses( mass_trees[0] );
    mass_tree_normalize_masses( mass_trees[1] );

    std::cout << "normalized KRD: " << earth_movers_distance( mass_trees[0], mass_trees[1] ) << std::endl;

    std::cout   << "normalized KRD, normalized by tree length: "
                << earth_movers_distance( mass_trees[0], mass_trees[1] ) / length( mass_trees[0] ) << std::endl;

    // ignore branch lengths: a more classical hop-based EMD
    mass_tree_transform_to_unit_branch_lengths( mass_trees[0] );
    mass_tree_transform_to_unit_branch_lengths( mass_trees[1] );

    std::cout << "normalized KRD, unit branch lengths: " << earth_movers_distance( mass_trees[0], mass_trees[1] ) << std::endl;

    std::cout << "total per edge error: " << calc_per_edge_error( scrapp_tree, true_tree ) << std::endl;

    return 0;
}
