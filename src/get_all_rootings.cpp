#include "genesis/genesis.hpp"

#include <string>
#include <sstream>
#include <vector>

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

void print_rootings(Tree const& tree,
										std::string const& out_dir)
{
	// for all edges
	for ( auto const& edge : tree.edges() ) {
		auto const id = edge.index();
		// get a rooted copy of the tree at current edge
		Tree cur_tree( tree );
		make_rooted( cur_tree, cur_tree.edge_at(id) );

		// print to file named after edge index
		std::stringstream ss;
		ss << out_dir << "edge_" << id << ".newick";
		// debug print
		// LOG_DBG << ss.str();
		// LOG_DBG << PrinterCompact().print( cur_tree );
		NewickWriter().to_file( cur_tree, ss.str() );
	}
}

int main(int argc, char* argv[]) {
  // Activate logging.
  utils::Logging::log_to_stdout();
  utils::Logging::details.date = true;
  utils::Logging::details.time = true;
  LOG_INFO << "Started";

	// Check if the command line contains the right number of arguments.
  if( argc != 3 ) {
      LOG_INFO << "Usage: " << argv[0] << " <tree_file> <output_dir>";
      return 1;
  }

	std::string ref_tree( argv[1] );
	Tree tree = NewickReader().read( from_file( ref_tree ) );

	std::string out_dir( argv[2] );
	if ( not is_dir( out_dir ) ) {
		throw std::runtime_error{std::string("output_dir doesn't exist or isnt a directory: ") + out_dir};
	}

	out_dir = dir_normalize_path( out_dir );

	// print the desired taxa lists to file, one per line
	print_rootings( tree, out_dir );

	LOG_INFO << "Finished";

	return 0;
}

