#include "genesis/genesis.hpp"

#include <string>
#include <sstream>
#include <vector>

using namespace genesis;
using namespace genesis::tree;
using namespace genesis::utils;

void print_rootings(Tree const& tree, std::string const& out_dir)
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
    CommonTreeNewickWriter().to_file( cur_tree, ss.str() );
  }
}

/*
  version of find_node that uses starts_with instead of equals
 */
TreeNode* find_node_starts_with(
    Tree& tree,
    const std::string& name
) {
  for( auto& node : tree.nodes() ) {
    if( starts_with( node.data<CommonNodeData>().name, name ) ) {
      return &node;
    }
  }

  return nullptr;
}

int main(int argc, char* argv[]) {
  // Activate logging.
  utils::Logging::log_to_stdout();
  utils::Logging::details.date = true;
  utils::Logging::details.time = true;
  LOG_INFO << "Started";

	// Check if the command line contains the right number of arguments.
  if( argc != 4 ) {
      LOG_INFO << "Usage: " << argv[0] << " <tree_file> <output_dir> <mode>";
      return 1;
  }

	std::string ref_tree( argv[1] );
	Tree tree = CommonTreeNewickReader().read( from_file( ref_tree ) );

	std::string out_dir( argv[2] );
	if ( not is_dir( out_dir ) ) {
		throw std::invalid_argument{std::string("output_dir doesn't exist or isn't a directory: ") + out_dir};
	}
	out_dir = dir_normalize_path( out_dir );

  std::string mode(argv[3]);

  if ( equals_ci( mode, "all") ) {
    // print all possible rootings of the tree
    print_rootings( tree, out_dir );
  } else if ( equals_ci( mode, "outgroup") ) {
    // find the outgroup and root by it / remove it
    auto names = node_names( tree );
    auto outgroup_node_ptr = find_node_starts_with( tree, "__SCRAPP_OUTGROUP__" );

    if ( outgroup_node_ptr == nullptr ) {
      throw std::runtime_error{"Unable to find outgroup taxon in file: " + ref_tree};
    }

    auto& new_root = outgroup_node_ptr->link().outer().node();

    change_rooting( tree, new_root );
    delete_leaf_node( tree, *outgroup_node_ptr );

    CommonTreeNewickWriter().to_file( tree, out_dir + "edge_byoutgroup.newick" );
  } else {
    throw std::invalid_argument{"Invalid mode: '" + mode + "'. Choose either 'all' or 'outgroup'."};
  }

	LOG_INFO << "Finished";

	return 0;
}

