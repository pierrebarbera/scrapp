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
    if( argc != 4 ) {
        LOG_INFO << "Usage: " << argv[0] << " <otu_file> <alignment_file> <output_file> ";
        return 1;
    }

    // In out paths.
    auto otu_file   = std::string( argv[1] );
    auto aln_file   = std::string( argv[2] );
    auto out_file   = std::string( argv[3] );

    auto fasta_reader = FastaReader();
    fasta_reader.site_casing( FastaReader::SiteCasing::kToLower );
    fasta_reader.guess_abundances( true );

    // read in otu file
    SequenceSet otus;
    fasta_reader.from_file( otu_file, otus );

    // STREAM in the aligned queries
    auto aligned_it = FastaInputIterator( fasta_reader ).from_file( aln_file );

    // prepare the output iterator
    std::ofstream out_stream( out_file );
    auto writer = FastaWriter().abundance_notation( FastaWriter::AbundanceNotation::kUnderscore );

    // for each aligned query
    while( aligned_it ) {
        // look up if the name exists in the otu file
        auto otu_seq = find_sequence( otus, aligned_it->label() );

        // if it exists there
        if ( otu_seq ) {
            Sequence out_seq( *aligned_it );
            // copy the multiplicity over
            out_seq.abundance( otu_seq->abundance() );
            // write to output
            writer.write_sequence( out_seq, out_stream);
        }

        aligned_it.increment();
    }

    LOG_INFO << "Finished";
    return 0;
}
