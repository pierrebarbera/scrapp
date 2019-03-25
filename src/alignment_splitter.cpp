/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech, Pierre Barbera

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

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <vector>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::sequence;
using namespace genesis::tree;
using namespace genesis::utils;

void merge_duplicates(
    SequenceSet& set
) {
    // Find duplicates, remove them and update the abundance counts of the originals
    std::unordered_map< std::string, size_t > dup_map;
    size_t i = 0;
    while( i < set.size() ) {
        auto& seq = set[i].sites();

        if( dup_map.count( seq ) == 0 ) {

            // If it is a new sequence, init the map entry and move to the next sequence.
            dup_map[ seq ] = i;
            ++i;

        } else {
            auto& original_seq = set[ dup_map[ seq ] ];
            original_seq.abundance( original_seq.abundance() + 1 );
            set.remove( i );
            // no need to iterate as we have effectively moved
        }
    }
}

static SequenceSet read_any_seqfile(std::string const& file)
{
    SequenceSet out_set;

    auto reader = PhylipReader();
    reader.site_casing( PhylipReader::SiteCasing::kToLower );
    try {
        reader.read( from_file( file ), out_set );
    } catch(std::exception& e) {
        out_set.clear();
        reader.mode( PhylipReader::Mode::kInterleaved );
        try {
            reader.read( from_file( file ), out_set );
        } catch(std::exception& e) {
            out_set.clear();
            try {
                FastaReader()
                    .site_casing( FastaReader::SiteCasing::kToLower )
                    .read( from_file( file ), out_set );
            } catch(std::exception& e) {
                throw std::invalid_argument{"Cannot parse sequence file(s): Invalid file format? (only phylip and fasta allowed)"};
            }
        }
    }

    return out_set;
}

std::vector<std::string> get_most_distant_leaf_per_node(
    Matrix<double> const& pwd,
    PlacementTree const& tree
) {
    auto ret = std::vector<std::string>( pwd.rows(), "" );

    // Nothing to do.
    if( pwd.rows() == 0 ) {
        return ret;
    }

    // search each row (whose indices equal the node indices) for the distance maximum
    // and the leaf behind it, then get that taxons name
    for( size_t r = 1; r < pwd.rows(); ++r ) {
        size_t best_c = 0;
        double cur_max = 0.0;
        for( size_t c = 0; c < pwd.cols(); ++c ) {
            if ( pwd( r, c ) > cur_max ) {
                // only act if this update to the current best would be a leaf
                // (this prevents an inner node to be selected as most distant in cases
                // where the leaf edge has dist = 0.0)
                if ( is_leaf( tree.node_at( best_c ) ) ) {
                    cur_max = pwd( r, c );
                    best_c = c;
                }
            }
        }
        ret[ r ] = tree.node_at( best_c ).data<PlacementNodeData>().name;
    }

    return ret;
}

// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if( argc < 4 or argc > 8 ) {
        LOG_INFO << "Usage: " << argv[0] << " <jplace_file> <query_alignment> <output_dir> "
                 << "<min_weight> <min_num> <max_num> <reference_alignment (activates outgroup inclusion)>";
        return 1;
    }

    // In out paths.
    auto jplace_file = std::string( argv[1] );
    auto query_alignment    = std::string( argv[2] );
    auto output_dir  = utils::trim_right( std::string( argv[3] ), "/") + "/";

    // If the argument is given, parse the min weight threshold.
    double min_weight = 0.51;
    size_t min_num = 4;
    size_t max_num = 500;
    bool include_outgroup = false;

    if( argc >= 5 ) {
        min_weight = std::atof( argv[4] );
    }
    // and the minimum number
    if( argc >= 6 ) {
        min_num = std::atoi( argv[5] );
    }
    // and the maximum number: if above this, do otu clustering
    if( argc >= 7 ) {
        max_num = std::atoi( argv[6] );
    }
    std::string reference_alignment;
    if( argc >= 8 ) {
        reference_alignment = argv[7];

        include_outgroup = file_exists( reference_alignment );

        if ( include_outgroup ) {
            min_num += 1;
        } else {
            throw std::invalid_argument{std::string("File '") + reference_alignment + "' does not exist."};
        }
    }

    LOG_INFO << "Specified: min_weight\t= " << min_weight;
    LOG_INFO << "Specified: min_num\t= " << min_num;

    // Read in placements.
    LOG_INFO << "Reading placements from jplace file.";
    auto sample = JplaceReader().read( from_file( jplace_file ) );
    LOG_DBG << "Found " << sample.size() << " Pqueries";

    rectify_values( sample );
    if( ! has_correct_edge_nums( sample.tree() )) {
        LOG_ERR << "Jplace edge nums were incorrect! Correcting them...";
        reset_edge_nums( sample.tree() );
    }

    // Only use the most probable placement.
    LOG_INFO << "Filtering placements.";
    filter_n_max_weight_placements( sample );

    // If the min weight is used, filter accordingly.
    if( min_weight < 1.0 ) {
        filter_min_weight_threshold( sample, min_weight );
        auto removed = remove_empty_pqueries( sample );
        LOG_DBG << "Filtered out " << removed << " Pqueries due to min weight threshold.";
    }

    // Get the pqueries per edge.
    auto const pqrs_per_edge = pqueries_per_edge( sample, true );
    assert( pqrs_per_edge.size() == sample.tree().edge_count() );

    // Read in sequences. Currently, we read the full file, because our Phylip reader does not yet
    // support to iterate a file (due to the stupid Phylip interleaved format). As we don't want to
    // maintain two implementations (Phylip with whole file and Fasta iteratievely), we simply
    // always read in the whole file. Might be worth improving in the future.
    SequenceSet seqs = read_any_seqfile( query_alignment );

    // Input check.
    if( ! is_alignment( seqs )) {
        LOG_ERR << "Alignment file contains sequences of different length, i.e., it is not an alignment.";
        return 1;
    }

    // dereplicate
    LOG_INFO << "Dereplicating sequences.";
    merge_duplicates( seqs );

    // after dereplication we add the ref seqs, if specified
    if ( include_outgroup ) {
        auto ref_seqs = read_any_seqfile( reference_alignment );
        for ( auto& s : ref_seqs ) {
            seqs.add( s );
        }
    }

    // Build map from seq label to id in the set.
    bool alignment_had_adundance = false;
    LOG_INFO << "Processing alignment.";
    auto seq_label_map = std::unordered_map< std::string, size_t >();
    for( size_t i = 0; i < seqs.size(); ++i ) {

        // Get label and check whether it is unique.
        auto const& label = seqs[i].label();
        if( seq_label_map.count( label ) > 0 ) {
            LOG_WARN << "Sequence '" << label << "' is not unique.";
        }

        if ( not alignment_had_adundance and seqs[i].abundance() > 1 ) {
            alignment_had_adundance = true;
        }

        seq_label_map[ label ] = i;
    }
    if ( alignment_had_adundance ) {
        LOG_INFO << "Detected abundance information in the original alignment, ignoring multiplicity count from the pqueries!";
    }

    // precompute node to node distance matrix to determine most distant taxon for each node
    std::vector<std::string> most_distant_taxon;
    if ( include_outgroup ) {
        auto pairwise_node_dists = node_branch_length_distance_matrix( sample.tree() );
        most_distant_taxon = get_most_distant_leaf_per_node( pairwise_node_dists, sample.tree() );
    }

    // Fill sequence sets for each edge.
    LOG_INFO << "Putting sequences on each branch.";
    auto edge_seqs = std::vector<SequenceSet>( sample.tree().edge_count() );
    auto finished_labels = std::unordered_map< std::string, size_t >();
    for( size_t edge_index = 0; edge_index < pqrs_per_edge.size(); ++edge_index ) {
        auto const& edge_pqueries = pqrs_per_edge[ edge_index ];

        // Get all names of the pqueries of the current edge.
        for( auto pqry_ptr : edge_pqueries ) {
            if( pqry_ptr->name_size() == 0 ) {
                LOG_WARN << "Pquery without name. Cannot extract sequence for this.";
            }
            if( pqry_ptr->name_size() > 1 ) {
                LOG_WARN << "Pquery with multiple names. This warning is just in case.";
            }

            for( auto const& pqry_name : pqry_ptr->names() ) {
                if( seq_label_map.count( pqry_name.name ) == 0 ) {
                    LOG_WARN << "No sequence found for Pquery '" << pqry_name.name << "'.";
                    continue;
                }

                // Add the sequence of the current name to the seq set of this edge.
                auto seq_id = seq_label_map[ pqry_name.name ];

                // if the abundance info wasnt in the original msa, and the seq hasnt had its abundance adjusted
                if ( not alignment_had_adundance
                    and seqs[ seq_id ].abundance() <= 1) {
                    // update abundance from jplace
                    seqs[ seq_id ].abundance( pqry_name.multiplicity );
                }

                edge_seqs[ edge_index ].add( seqs[ seq_id ] );
                ++finished_labels[ pqry_name.name ];
            }
        }

        // if outgroup inclusion was selected, add that taxon to the set (with a little renaming for clarity)
        if ( include_outgroup and not most_distant_taxon.empty() ) {
            auto const& edge = sample.tree().edge_at( edge_index );
            auto const primary_node_id = edge.primary_node().index();

            auto seq_id = seq_label_map[ most_distant_taxon[ primary_node_id ] ];

            Sequence seq_cpy = seqs[ seq_id ];
            seq_cpy.label( "__SCRAPP_OUTGROUP__" + seq_cpy.label() );

            edge_seqs[ edge_index ].add( seq_cpy );
        }
    }

    // Sanity checks.
    for( auto const& label_count : finished_labels ) {
        if( label_count.second > 1 ) {
            LOG_WARN << "Encountered label more than once: '" << label_count.first << "'.";
        }
        seq_label_map.erase( label_count.first );
    }
    if( seq_label_map.size() > 0 ) {
        LOG_WARN << "There were " << seq_label_map.size() << " sequences that did not occur "
                 << "in any Pquery. Those might be simply the original reference sequences. "
                 << "In case these numbers to match, better check if everything is correct.";
    }

    // Write result files.
    LOG_INFO << "Writing result files.";
    auto writer = FastaWriter();
    // ensure abundance info is written out
    writer.abundance_notation( FastaWriter::AbundanceNotation::kUnderscore );
    for( size_t edge_index = 0; edge_index < edge_seqs.size(); ++edge_index ) {

        // Check: Don't need to write empty sequence files.
        if( edge_seqs[ edge_index ].size() < min_num ) {
            continue;
        }

        // !! edge_num isnt the same ad the edge_index! former is as in jplace, latter is genesis internal
        auto const& edge = sample.tree().edge_at( edge_index );
        auto const edge_num = edge.data<PlacementEdgeData>().edge_num();

        LOG_INFO << "\tEdge: " << edge_num << "    \tSize: " << edge_seqs[ edge_index ].size();

        // create edge outdir
        auto edge_dir = output_dir + "edge_" + std::to_string( edge_num ) + "/";
        dir_create( edge_dir );

        // Write result
        writer.to_file( edge_seqs[ edge_index ], edge_dir + "aln.fasta" );

        if( edge_seqs[ edge_index ].size() > max_num ) {
            LOG_INFO << "\t\tExceeded max number of sequences, writing stripped file to signal clustering";
            remove_all_gaps( edge_seqs[ edge_index ] );
            writer.to_file( edge_seqs[ edge_index ], edge_dir + "stripped.fasta" );
        }
    }

    LOG_INFO << "Finished";
    return 0;
}
