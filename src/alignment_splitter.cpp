/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2019 Pierre Barbera and Lucas Czech

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

namespace here {

    class Range {
    public:
      Range (const size_t begin, const size_t end)
        : begin(begin), end(end) {assert(end >= begin);};
      Range() = default;
      ~Range () = default;

      operator bool() const { return end > 0;}

      size_t begin;
      size_t end;

      size_t span() { return end - begin;};

    private:

    };

    here::Range get_valid_range(const std::string& sequence)
    {
      size_t begin = 0;
      size_t end = sequence.length();

      assert(end);

      while(begin < end and sequence.c_str()[begin] == '-') {
        begin++;
      }

      while(end > begin and sequence.c_str()[end - 1u] == '-') {
        end--;
      }

      return Range(begin, end);
    }

    int overlap( const here::Range& lhs, const here::Range& rhs )
    {
        auto begin = std::max(lhs.begin, rhs.begin);
        auto end = std::min(lhs.end, rhs.end);
        return end - begin;
    }

}

void merge_duplicates( SequenceSet& set )
{
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
    for( size_t r = 0; r < pwd.rows(); ++r ) {
        size_t best_c = 0;
        double cur_max = 0.0;
        for( size_t c = 0; c < pwd.cols(); ++c ) {
            if ( pwd( r, c ) > cur_max ) {
                // only act if this update to the current best would be a leaf
                // (this prevents an inner node to be selected as most distant in cases
                // where the leaf edge has dist = 0.0)
                if ( is_leaf( tree.node_at( c ) ) ) {
                    cur_max = pwd( r, c );
                    best_c = c;
                }
            }
        }

        assert ( is_leaf( tree.node_at( best_c ) ) );

        ret[ r ] = tree.node_at( best_c ).data<PlacementNodeData>().name;
    }

    return ret;
}

void overlap_check( SequenceSet const& set, int const min_overlap, double const overlap_ratio_thresh )
{
    std::vector<here::Range> ranges;
    ranges.reserve( set.size() );

    for (size_t i = 0; i < set.size(); ++i) {
        auto& seq = set[i].sites();

        // get the valid range
        ranges.emplace_back( here::get_valid_range( seq ) );
    }

    std::sort( ranges.begin(), ranges.end(),
        []( const here::Range& lhs, const here::Range& rhs ){
            if ( lhs.begin < rhs.begin ) {
                return true;
            } else if ( lhs.begin == rhs.begin ) {
                // tiebreaker needed
                return lhs.end > rhs.end;
            } else {
                return false;
            }
        }
    );

    for (size_t i = 0; i < ranges.size() - 1u; ++i) {
        auto j = i + 1u;
        if ( here::overlap(ranges[i], ranges[j]) < min_overlap ) {
            throw std::runtime_error{"Reads do not sufficiently overlap!"};
        }
    }

    // pairwise overlap percentage check
    std::vector<double> ratios;
    for (size_t i = 0; i < ranges.size() - 1u; ++i) {
        for (size_t j = i + 1u; j < ranges.size(); ++j) {
            ratios.push_back( static_cast<double>(here::overlap(ranges[i], ranges[j]))
                / std::min(ranges[i].span(), ranges[j].span()) );
        }
    }

    double avg = std::accumulate(ratios.begin(), ratios.end(), 0.0) / ratios.size();

    if ( avg < overlap_ratio_thresh ) {
        throw std::runtime_error{"Average pairwise read overlap below threshold! Average: " + std::to_string(avg)};
    }

}

using Pquery_set = std::vector< Pquery const* >;
using Point = std::vector<double>;

class sortable_pquery
{
public:
    sortable_pquery( Pquery const* ptr ) : pq_ptr( ptr ) {};
    ~sortable_pquery() = default;
    Pquery const* pq_ptr = nullptr;
    size_t sort_key = 0;
    bool operator < ( sortable_pquery const& other ) const {
        return sort_key < other.sort_key;
    }
    bool operator==( sortable_pquery const& other ) const {
        return pq_ptr == other.pq_ptr;
    }

    operator Pquery const* () const {
        return pq_ptr;
    }
};

PqueryPlacement const& best_placement( Pquery const& pq ) {
    assert( pq.placement_size() != 0 );

    PqueryPlacement const * best = &pq.placement_at(0);
    for( auto const& p : pq.placements() ) {
        if( best->like_weight_ratio < p.like_weight_ratio ) {
            best = &p;
        }
    }

    return *best;
}

double distance( Point const& lhs, Point const& rhs, size_t const dimensions=2 )
{
    // Simple euclidean distance
    double sum = 0.0;
    for( size_t d = 0; d < dimensions; ++d ) {
        auto const diff = lhs[d] - rhs[d];
        sum += diff * diff;
    }
    return std::sqrt( sum );
}

Sequence const& get_sequence( SequenceSet const& haystack, std::string const& needle )
{
    auto iter = std::find_if(std::begin(haystack), std::end(haystack),
                [&needle](Sequence const& s){
                    return s.label() == needle;
    });

    if( iter == std::end(haystack) ) {
        throw std::invalid_argument{
            std::string("Sequence of name '") + needle +
            "' doesn't have a sequence in the input fasta file."
        };
    }

    return *iter;
}

std::string const& name( Pquery const* pq_ptr )
{
    return pq_ptr->name_at( 0 ).name;
}

size_t count_informative_sites( Sequence const& seq )
{
    auto relevant_chars = nucleic_acid_codes_plain() + nucleic_acid_codes_degenerated();

    // make a counter for each value of char
    std::vector<size_t> counters( std::numeric_limits< char >::max() + 1, 0 );

    for( auto const& c : seq ) {
        counters[ c ]++;
    }

    // sum up only those chars we are interested int
    size_t counter = 0;
    for( auto const& c : relevant_chars ) {
        counter += counters[ c ];
    }

    return counter;
}


Pquery_set placement_pseudokmedoids( Pquery_set const& set, SequenceSet const& msa, size_t const k, size_t const x )
{
    assert( k >= 2 );

    // build a vector that EuclidianKmeans can work with
    // aka translating placements into the placement space! hell yes

    std::vector<Point> placement_coordinates;
    placement_coordinates.reserve( set.size() );

    for( auto pq : set ) {
        // get best placement
        auto const& p = best_placement( *pq );

        // add this pquery's placement coordinates to thed dataset
        placement_coordinates.push_back( { p.pendant_length, p.proximal_length } );
    }

    // do the kmeans
    auto clustering = EuclideanKmeans( 2 );
    clustering.run( placement_coordinates, k );
    auto centroids = clustering.centroids();
    auto const& point_id_to_cluster = clustering.assignments();

    // make a set of clusters
    assert( point_id_to_cluster.size() == set.size() );
    std::vector<std::vector< sortable_pquery >> clusters( k );
    for( size_t point_id = 0; point_id < point_id_to_cluster.size(); ++point_id ) {
        auto const cluster_id = point_id_to_cluster[ point_id ];
        clusters[ cluster_id ].push_back( set[ point_id ] );
    }

    // find the closest medoid to every centroid
    // assert( cluster_reps.size() == centroids.size() );
    assert( k == centroids.size() );

    // pick the top x sequences per cluster as the cluster representatives
    for( auto& cluster : clusters ) {
        // ensure x is strictly smaller. in all other cases, keep all representatives
        if(  x >= cluster.size() ) {
            continue;
        }

        for( auto& spq : cluster ) {
            // find that sequence in the query msa
            auto const& seq = get_sequence( msa, name( spq.pq_ptr ) );

            // set this pqueries sort key: the number of informative sites
            spq.sort_key = count_informative_sites( seq );
        }

        // sort the datapoints of this cluster (ascending!)
        std::sort( std::begin(cluster), std::end(cluster) );

        // remove all but the best x elements from the container
        auto remove_iter = std::end(cluster);
        std::advance( remove_iter, -x );
        cluster.erase( std::begin(cluster), remove_iter );
    }

    Pquery_set cluster_reps;
    cluster_reps.reserve( k * x );
    for( auto const& cluster : clusters ) {
        cluster_reps.insert( std::end(cluster_reps), std::begin(cluster), std::end(cluster) );
    }


    // check that cluster_reps are unique
    for( auto pq_ptr : cluster_reps ) {
        auto num = std::count( std::begin(cluster_reps), std::end(cluster_reps), pq_ptr );
        if( num > 1 ) {
            std::cout << "ERR: Representative '" << name( pq_ptr ) << "' occured " << num << " times!" << std::endl;
            throw std::runtime_error{"Abort, see above"};
        }
    }

    return cluster_reps;
}
// =================================================================================================
//     Main
// =================================================================================================

int main( int argc, char** argv )
{

    utils::Options::get().allow_file_overwriting( true );
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

    merge_duplicates( sample );

    // Get the pqueries per edge.
    auto pqrs_per_edge = pqueries_per_edge( sample, true );
    assert( pqrs_per_edge.size() == sample.tree().edge_count() );

    // Read in sequences. Currently, we read the full file, because our Phylip reader does not yet
    // support to iterate a file (due to the stupid Phylip interleaved format). As we don't want to
    // maintain two implementations (Phylip with whole file and Fasta iteratievely), we simply
    // always read in the whole file. Might be worth improving in the future.
    SequenceSet seqs = read_any_seqfile( query_alignment );

    // Input check.
    if( not is_alignment( seqs )) {
        LOG_ERR << "Alignment file contains sequences of different length, i.e., it is not an alignment.";
        return 1;
    }

    // dereplicate
    // LOG_INFO << "Dereplicating sequences.";
    // merge_duplicates( seqs );

    // after dereplication we add the ref seqs, if specified
    SequenceSet ref_seqs;
    if ( include_outgroup ) {
        ref_seqs = read_any_seqfile( reference_alignment );
        if ( ref_seqs.empty() ) {
            LOG_ERR << "Reference alignment file was empty.";
            return 1;
        } else if( not is_alignment( ref_seqs )) {
            LOG_ERR << "Reference alignment file contains sequences of different length, i.e., it is not an alignment.";
            return 1;
        }  else if( seqs.at( 0 ).length() != ref_seqs.at( 0 ).length() ) {
            LOG_ERR << "Reference alignment and Query alignment have differing widths.";
            return 1;
        }
        LOG_INFO << "Outgroup inclusion was selected, reference MSA found.";
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

        if( most_distant_taxon.empty()) {
            throw std::runtime_error{"get_most_distant_leaf_per_node failed?"};
        }
    }

    // Fill sequence sets for each edge.
    LOG_INFO << "Putting sequences on each branch.";
    auto edge_seqs = std::vector<SequenceSet>( sample.tree().edge_count() );
    auto finished_labels = std::unordered_map< std::string, size_t >();
    for( size_t edge_index = 0; edge_index < pqrs_per_edge.size(); ++edge_index ) {
        auto& edge_pqueries = pqrs_per_edge[ edge_index ];

        // IF there are too many pqueries on this edge, perform clustering
        if( edge_pqueries.size() > max_num ) {
            size_t const k = max_num/10;
            size_t const x = 10;
            LOG_INFO << "Number of queries above threshold! clustering. k = " << k <<", x = " << x;
            edge_pqueries = placement_pseudokmedoids( edge_pqueries, seqs, k, x );
        }

        // Get all names of the pqueries of the current edge.
        for( auto pqry_ptr : edge_pqueries ) {
            if( pqry_ptr->name_size() == 0 ) {
                LOG_WARN << "Pquery without name. Cannot extract sequence for this.";
            }
            if( pqry_ptr->name_size() > 1 ) {
                LOG_WARN << "Pquery with multiple names. This warning is just in case.";
            }

            // Find the sequence associated with the first label of the pquery (we only want one representative)
            auto seq_id = seq_label_map[ pqry_ptr->name_at(0).name ];

            for( auto const& pqry_name : pqry_ptr->names() ) {
                if( seq_label_map.count( pqry_name.name ) == 0 ) {
                    LOG_WARN << "No sequence found for Pquery '" << pqry_name.name << "'.";
                    continue;
                }

                // update abundance accounting for current extra sequence label
                auto abund = seqs[ seq_id ].abundance();
                seqs[ seq_id ].abundance( abund + pqry_name.multiplicity );

                ++finished_labels[ pqry_name.name ];
            }
            edge_seqs[ edge_index ].add( seqs[ seq_id ] );
        }

        // if outgroup inclusion was selected, add that taxon to the set (with a little renaming for clarity)
        if ( include_outgroup and not most_distant_taxon.empty() ) {
            auto const& edge = sample.tree().edge_at( edge_index );
            auto const primary_node_id = edge.primary_node().index();

            Sequence seq_cpy = get_sequence( ref_seqs, most_distant_taxon[ primary_node_id ] );
            LOG_INFO << "Including sequence as outgroup: '" << seq_cpy.label() << "'";
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

        auto& seqs = edge_seqs[ edge_index ];

        // Check: Don't need to write empty sequence files.
        if( seqs.size() < min_num ) {
            continue;
        }

        // !! edge_num isnt the same ad the edge_index! former is as in jplace, latter is genesis internal
        auto const& edge = sample.tree().edge_at( edge_index );
        auto const edge_num = edge.data<PlacementEdgeData>().edge_num();

        LOG_INFO << "\tEdge: " << edge_num << "    \tSize: " << edge_seqs[ edge_index ].size();

        overlap_check( seqs, 50, 0.5 );

        // create edge outdir
        auto edge_dir = output_dir + "edge_" + std::to_string( edge_num ) + "/";
        dir_create( edge_dir );

        // Write result
        writer.to_file( seqs, edge_dir + "aln.fasta" );

    }

    LOG_INFO << "Finished";
    return 0;
}
