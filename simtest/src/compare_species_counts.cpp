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

using per_edge_func = std::function<double(double,double)>;

std::vector<double> calc_per_edge_errors( Tree const& lhs, Tree const& rhs, per_edge_func& func, bool allow_zero=true )
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

        if ( rhs_sc != 0.0 or allow_zero ) {
            auto diff = func(lhs_sc, rhs_sc);

            per_edge_error.emplace_back( diff );
        }
    }

    // Now we need to be done with both trees, otherwise we have a problem.
    if( lhs_it != lhs_end || rhs_it != rhs_end ) {
        throw std::invalid_argument( "Incompatible MassTrees." );
    }

    // sort so that we may calc the median later
    std::sort( std::begin(per_edge_error), std::end(per_edge_error) );

    return per_edge_error;
}

double total_counts( Tree const& tree )
{
    double total = 0.0;
    for( auto& edge : tree.edges() ) {
        total += edge.data<ScrappEdgeData>().species_count;
    }
    return total;
}

double total_counts_exp( Tree const& tree )
{
    double total = 0.0;
    for( auto& edge : tree.edges() ) {
        total += std::exp( edge.data<ScrappEdgeData>().species_count) ;
    }
    return total;
}

void normalize_species_counts( Tree& tree )
{
    auto total = total_counts( tree );

    for( auto& edge : tree.edges() ) {
        edge.data<ScrappEdgeData>().species_count = ( edge.data<ScrappEdgeData>().species_count ) / total;
    }
}

void modify_species_counts( Tree& tree )
{
    for( auto& edge : tree.edges() ) {
        edge.data<ScrappEdgeData>().species_count = edge.data<ScrappEdgeData>().species_count * 1.0;
    }
}

// =================================================================================================
//     main
// =================================================================================================

int main( int argc, char** argv )
{

    // Check if the command line contains the right number of arguments.
    if( argc != 4 ) {
      LOG_INFO << "Usage: " << argv[0] << " <scrapp_result (newick)> <true_annotation (newick)> <csv_mode>";
      return 1;
    }

    std::string scrapp_tree_file( argv[1] );
    std::string true_tree_file( argv[2] );
    bool const csv_mode = std::atoi( argv[3] );
    char const sep = ',';

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

    auto const krd = earth_movers_distance( mass_trees[0], mass_trees[1] );

    if ( csv_mode ) {
        std::cout << krd << sep;
    } else {
        std::cout << "pure KRD: " << krd << std::endl;
    }

    // mass normalization
    mass_tree_normalize_masses( mass_trees[0] );
    mass_tree_normalize_masses( mass_trees[1] );

    auto const norm_krd = earth_movers_distance( mass_trees[0], mass_trees[1] );
    auto const norm_norm_krd = norm_krd / length( mass_trees[0] );

    if ( csv_mode ) {
        std::cout   << norm_krd << sep;
        std::cout   << norm_norm_krd << sep;
    } else {
        std::cout   << "normalized KRD: " << norm_krd << std::endl;
        std::cout   << "normalized KRD, normalized by tree length: "
                    << norm_norm_krd << std::endl;
    }

    // ignore branch lengths: a more classical hop-based EMD
    mass_tree_transform_to_unit_branch_lengths( mass_trees[0] );
    mass_tree_transform_to_unit_branch_lengths( mass_trees[1] );

    auto const norm_unit_krd = earth_movers_distance( mass_trees[0], mass_trees[1] );

    if ( csv_mode ) {
        std::cout << norm_unit_krd << sep;
    } else {
        std::cout << "normalized KRD, unit branch lengths: " << norm_unit_krd << std::endl;
    }

    modify_species_counts( scrapp_tree );

    per_edge_func abs_diff_func = [](double l, double r){return std::abs(l - r);};

    auto per_edge_errors = calc_per_edge_errors( scrapp_tree, true_tree, abs_diff_func);

    per_edge_func diff_func = [](double l, double r){return l - r;};
    auto per_edge_diff = calc_per_edge_errors( scrapp_tree, true_tree, diff_func);
    HistogramAccumulator histacc( per_edge_diff );
    std::ofstream hist("hist");
    hist << histacc.build_uniform_ranges_histogram( 10 );


    auto const total_true_count = total_counts( true_tree );
    auto const num_counts = per_edge_errors.size();

    auto error_sum = std::accumulate( std::begin(per_edge_errors), std::end(per_edge_errors), 0.0);
    auto minmax = minimum_maximum( std::begin(per_edge_errors), std::end(per_edge_errors) );
    auto meanstddev = mean_stddev( std::begin(per_edge_errors), std::end(per_edge_errors) );
    auto med = median( std::begin(per_edge_errors), std::end(per_edge_errors) );

    if ( csv_mode ) {
        std::cout << error_sum << sep << meanstddev.mean << sep << med << sep;
    } else {
        std::cout << "~~~ per edge errors ~~~" << std::endl;
        std::cout << "\ttottru:\t" << total_counts( true_tree ) << std::endl;
        std::cout << "\ttotinf:\t" << total_counts( scrapp_tree ) << std::endl;
        std::cout << "\treldev:\t" << meanstddev.mean / (total_true_count / num_counts) << std::endl;
        std::cout << "\tsum:\t" << error_sum << std::endl;
        std::cout << "\tmean:\t" << meanstddev.mean << std::endl;
        std::cout << "\tstddev:\t" << meanstddev.stddev << std::endl;
        std::cout << "\tmedian:\t" << med << std::endl;
        std::cout << "\tmin:\t" << minmax.min << std::endl;
        std::cout << "\tmax:\t" << minmax.max << std::endl;
    }

    per_edge_func rel_diff_func = [](double l, double r){return std::abs((l - r)/r);};
    per_edge_errors = calc_per_edge_errors( scrapp_tree, true_tree, rel_diff_func, false);


    error_sum = std::accumulate( std::begin(per_edge_errors), std::end(per_edge_errors), 0.0);
    minmax = minimum_maximum( std::begin(per_edge_errors), std::end(per_edge_errors) );
    meanstddev = mean_stddev( std::begin(per_edge_errors), std::end(per_edge_errors) );
    med = median( std::begin(per_edge_errors), std::end(per_edge_errors) );

    if ( csv_mode ) {
        std::cout << error_sum << sep << meanstddev.mean << sep << med << sep;
    } else {
        std::cout << "~~~ using relative error ~~~" << std::endl;
        std::cout << "\ttottru:\t" << total_counts( true_tree ) << std::endl;
        std::cout << "\ttotinf:\t" << total_counts( scrapp_tree ) << std::endl;
        std::cout << "\treldev:\t" << meanstddev.mean / (1.0 / num_counts) << std::endl;
        std::cout << "\tsum:\t" << error_sum << std::endl;
        std::cout << "\tmean:\t" << meanstddev.mean << std::endl;
        std::cout << "\tstddev:\t" << meanstddev.stddev << std::endl;
        std::cout << "\tmedian:\t" << med << std::endl;
        std::cout << "\tmin:\t" << minmax.min << std::endl;
        std::cout << "\tmax:\t" << minmax.max << std::endl;
    }


    // normalize the species counts per tree
    normalize_species_counts( scrapp_tree );
    normalize_species_counts( true_tree );

    per_edge_errors = calc_per_edge_errors( scrapp_tree, true_tree, abs_diff_func );


    error_sum = std::accumulate( std::begin(per_edge_errors), std::end(per_edge_errors), 0.0);
    minmax = minimum_maximum( std::begin(per_edge_errors), std::end(per_edge_errors) );
    meanstddev = mean_stddev( std::begin(per_edge_errors), std::end(per_edge_errors) );
    med = median( std::begin(per_edge_errors), std::end(per_edge_errors) );

    if ( csv_mode ) {
        std::cout << error_sum << sep << meanstddev.mean << sep << med;
    } else {
        std::cout << "~~~ after normalization ~~~" << std::endl;
        std::cout << "\ttottru:\t" << total_counts( true_tree ) << std::endl;
        std::cout << "\ttotinf:\t" << total_counts( scrapp_tree ) << std::endl;
        std::cout << "\treldev:\t" << meanstddev.mean / (1.0 / num_counts) << std::endl;
        std::cout << "\tsum:\t" << error_sum << std::endl;
        std::cout << "\tmean:\t" << meanstddev.mean << std::endl;
        std::cout << "\tstddev:\t" << meanstddev.stddev << std::endl;
        std::cout << "\tmedian:\t" << med << std::endl;
        std::cout << "\tmin:\t" << minmax.min << std::endl;
        std::cout << "\tmax:\t" << minmax.max << std::endl;
    }

    if ( csv_mode ) {
        std::cout << std::endl;
    }

    return 0;
}
