#include <iostream>
#include <vector>
#include <map>
#include <assert.h>
#include <stdlib.h> // to generate random number
#include <ctime>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <fstream>
#include <chrono>
#include <string>

#pragma once

class gen_topology
{

private:
	int kr; // number of rows in baseline mesh topology
	int kc; // number of column in baseline mesh topology
	int links_removed;  // number of uni-directional links to be removed
	int num_topology;  // number of topologies to be generated
	std::vector< std::vector<int> > mesh;
	// this contains link between 'src'--'dest' pairs
	std::vector< std::vector<int> > connectvty_mtrx;
	// std::map<int, int> connectvty_map;
	// this is new topology
	// This structure will check if the given 'coonectvty_matrix' is connected
	std::vector<bool> visited;
public:
	std::vector< std::vector<int> > new_connectvty_mtrx; // this is the new topology
	std::vector< std::vector<int> > all_pairs_shortest_path;
    std::vector< std::vector<std::vector<std::vector<int> > > > all_pairs_all_paths;
	std::vector< std::vector<char> > up_down_mtrx;
	gen_topology()
		: kr(0)
		, kc(0)
		, links_removed(0)
		, num_topology(0)
	{ init(); }

	gen_topology(int kr_, int kc_, int links_removed_, int num_topology_)
			: kr(kr_)
			, kc(kc_)
			, links_removed(links_removed_)
			, num_topology(num_topology_)
			{ init(); }

	~gen_topology() {
		delete this;
	}
	void print_matrix(std::vector< std::vector<int> >& matrx);
	void print_char_matrix(std::vector< std::vector<char> >& matrx);
	bool check_connected();
	void populate_all_pair_all_paths();
	void findAllPaths(int src, int dst);
	bool findAllPathsUtil(int, int, std::vector<bool>&,
	                        std::vector<int>&, int&, int&);
	void DFS_(int start_node);
	void prnt_state();
	void init();
	void generate_topology();
	bool remove_link(int src_node);
	void populate_up_dwn_mtrx();
	inline int get_num_topology() {
		return num_topology;
	}
	inline int get_kr() {
		return kr;
	}
	inline int get_kc() {
		return kc;
	}
	inline int get_links_removed() {
		return links_removed;
	}
};
	void usage();