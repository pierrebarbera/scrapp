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

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::sequence;
using namespace genesis::tree;
using namespace genesis::utils;

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
    if( argc < 4 or argc > 6 ) {
        LOG_INFO << "Usage: " << argv[0] << " <jplace_file> <alignment_file> <output_dir> "
                 << "<min_weight> <min_num>";
        return 1;
    }

    // In out paths.
    auto jplace_file = std::string( argv[1] );
    auto aln_file    = std::string( argv[2] );
    auto output_dir  = utils::trim_right( std::string( argv[3] ), "/") + "/";

    // If the argument is given, parse the min weight threshold.
    double min_weight = 0.51;
    size_t min_num = 4;
    if( argc >= 5 ) {
        min_weight = std::atof( argv[4] );
    }
    if( argc == 6 ) {
        min_num = std::atoi( argv[5] );
    }

    LOG_INFO << "Specified: min_weight\t= " << min_weight;
    LOG_INFO << "Specified: min_num\t= " << min_num;

    // Read in placements.
    LOG_INFO << "Reading placements from jplace file.";
    auto sample = JplaceReader().from_file( jplace_file );
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
    SequenceSet seqs;
    try{
        LOG_INFO << "Reading alignment as Phylip file.";
        auto phylip_reader = PhylipReader();
        phylip_reader.site_casing( PhylipReader::SiteCasing::kToLower );
        phylip_reader.from_file( aln_file, seqs );
    } catch( ... ) {
        LOG_INFO << "Phylip failed, trying Fasta now.";
        LOG_INFO << "Reading alignment as Fasta file.";
        seqs.clear();
        auto fasta_reader = FastaReader();
        fasta_reader.site_casing( FastaReader::SiteCasing::kToLower );
        fasta_reader.from_file( aln_file, seqs );
    }

    // Input check.
    if( ! is_alignment( seqs )) {
        LOG_ERR << "Alignment file contains sequences of different length, i.e., it is not an alignment.";
        return 1;
    }

    // dereplicate
    LOG_INFO << "Dereplicating sequences.";
    merge_duplicate_sequences( seqs );

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

    LOG_INFO << "De-aligning sequences.";

    for( auto& seq_set : edge_seqs) {
        // merge_duplicate_sequences( seq_set, MergeDuplicateSequencesCountPolicy::kAppendToLabel );
        remove_all_gaps( seq_set );
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

        LOG_INFO << "\tEdge: " << edge_index << "    \tSize: " << edge_seqs[ edge_index ].size();

        // create edge outdir
        auto edge_dir = output_dir + "edge_" + std::to_string( edge_index ) + "/";
        dir_create(edge_dir);

        // Write to Phylip.
        std::string output_file = edge_dir + "aln.fasta";
        writer.to_file( edge_seqs[ edge_index ], output_file );
    }

    LOG_INFO << "Finished";
    return 0;
}
