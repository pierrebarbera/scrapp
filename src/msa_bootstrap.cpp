#include "genesis/genesis.hpp"

#include <string>
#include <random>
#include <sstream>


using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;
using namespace genesis::sequence;

int main(int argc, char* argv[]) {
    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.date = true;
    utils::Logging::details.time = true;
    LOG_INFO << "Started";

	// Check if the command line contains the right number of arguments.
    if( argc != 3 ) {
        LOG_INFO << "Usage: " << argv[0] << " <msa_file> <output_dir> <num_reps>";
        return 1;
    }

	std::string aln_file( argv[1] );

	std::string out_dir( argv[2] );
	if ( not is_dir( out_dir ) ) {
		throw std::runtime_error{std::string("output_dir doesn't exist or isn't a directory: ") + out_dir};
	}
	out_dir = dir_normalize_path( out_dir );

    int const num_reps = std::atoi( argv[3] );
    if( num_reps < 1 ) {
        throw std::runtime_error{
            std::string("Invalid number of replicates: ") + std::to_string( num_reps )
        };
    }

    // read MSA
    SequenceSet seqs;
    try{
        LOG_INFO << "Reading alignment as Phylip file.";
        auto phylip_reader = PhylipReader();
        phylip_reader.site_casing( PhylipReader::SiteCasing::kToLower );
        phylip_reader.read( from_file( aln_file ), seqs );
    } catch( ... ) {
        LOG_INFO << "Phylip failed, trying Fasta now.";
        LOG_INFO << "Reading alignment as Fasta file.";
        seqs.clear();
        auto fasta_reader = FastaReader();
        fasta_reader.site_casing( FastaReader::SiteCasing::kToLower );
        fasta_reader.read( from_file( aln_file ), seqs );
    }

    if ( seqs.empty() ) {
        throw std::runtime_error{std::string("File is empty: ") + aln_file};
    }

    if ( not is_alignment( seqs ) ) {
        throw std::runtime_error{std::string("File isn't a MSA: ") + aln_file};
    }


    size_t const width = seqs.at( 0 ).size();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, width - 1u);

    // generate desired number of MSA replicates
    for (size_t replicate = 0; replicate < static_cast<size_t>( num_reps ); ++replicate) {
        std::ofstream cur_out_file( out_dir + "replicate_" + std::to_string( replicate ) + ".fasta" );
        FastaOutputIterator out( cur_out_file );

        // by first generating a random sequence of columns
        std::vector<size_t> columns( width );
        for (size_t i = 0; i < width; ++i) {
            columns[i] = dis(gen);
        }

        // then generating the sequences on write
        for ( auto const& row : seqs ) {
            // make a temp sequence for output
            Sequence rand_seq( row.label(), std::string( row.size(), '0' ), row.abundance() );
            for ( size_t i = 0; i < row.size(); ++i ) {
                rand_seq.sites()[ i ] = row[ columns[ i ] ];
            }
            out = rand_seq;
        }
    }

	LOG_INFO << "Finished";

	return 0;
}

