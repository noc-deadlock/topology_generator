#include "gen_topology.hh"

using namespace std;

// Globals
/******************/
int num_row = 3;
int num_col = 3;
int link_removed = 0;
int num_topology = 1;
/*****************/
void usage (void);

/******************************
 * Helper function to print 2D int Matrix
 * ****************************/
void
gen_topology::print_matrix(std::vector< std::vector<int> >& matrx) {

	for(unsigned row = 0; row < matrx.size(); ++row) {
		cout << "Row ["<<row<<"]:\t";
		for(unsigned col = 0; col < matrx.size(); ++col) {
			cout << matrx[row][col] << " ";
		}
		cout << endl;
	}
}

/******************************
 * Helper function to print 2D char Matrix
 * ****************************/
void
gen_topology::print_char_matrix(std::vector<std::vector<char> > &matrx) {

	for(unsigned i = 0; i < matrx.size(); ++i) {
		for(unsigned j = 0; j < matrx.size(); ++j) {
			cout << matrx[i][j] << " ";
		}
		cout << endl;
	}
}
/*******************************
Does recursive DFS_ search to find the
 topology to be connected.
********************************/
// all member definition goes here
void
gen_topology::DFS_(int v) {
	// mark the current node as visited and print it
	visited[v] = true;
	std::vector<int> row_ = new_connectvty_mtrx[v];
	for(unsigned i = 0; i < row_.size(); ++i) {
		if(row_[i] == 1) {
			if(visited[i] == false) {
				DFS_(i);
			}
		}
	}
}

/*******************************
Top-level function to check if the
 resulting topology is connected or not.
********************************/
/*
Return false if disconnected... returns true otherwise..
*/
bool
gen_topology::check_connected() {
	// 1. Start by clearing visited vector
	for(unsigned i = 0; i < visited.size(); ++i) {
		visited[i] = false;
	}
	// 2. Call the DFS_(0); yu can always start with node 0
	// doesn't matter
	DFS_(0); // this will populate the 'visited' vector
	// 3. is any element of visited-vector is false.. return
	// false otherwise true
	for(unsigned i = 0; i < visited.size(); ++i) {
		if(visited[i] == false) {
			return false;
		}
	}

	return true;
}

void
gen_topology::prnt_state() {

	cout << "kr: " << kr << endl;
	cout << "kc: " << kc << endl;
	cout << "links_removed: " << links_removed << endl;
	cout << "num_topology: " << num_topology << endl;
	cout << "-------Mesh--------" << endl;

	for(int row = kr-1; row >= 0; --row) {
		cout << "Row ["<<row<<"]:\t";
		for(int col = 0; col < kc; ++col) {
			cout << mesh[row][col] << " ";
		}
		cout << endl;
	}

	cout << "-------connectvty_mtrx--------" << endl;
	for(int row = 0; row < (kr*kc); ++row) {
		cout << "Row ["<<row<<"]:\t";
		for(int col = 0; col < (kr*kc); ++col) {
			cout << connectvty_mtrx[row][col] << " ";
		}
		cout << endl;
	}

	return;
};

/*******************************
This function initializes the connectivity
 matrix of baseline Mesh topology
********************************/
void
gen_topology::init() {
	// initialize the mesh..
	// row => y-coordinate
	// col => x-coordinate
	// to access any coordinate: mesh[y][x] => mesh[row][col]
	mesh.resize(kr);
	for(int i = 0; i < kr; ++i) {
		mesh[i].resize(kc);
	}
	// Mesh contains the node ids..
	int id_ = 0;
	for(int row = 0; row < kr; ++row) {
		for(int col = 0; col < kc; ++col) {
			mesh[row][col] = id_++;
		}
	}
	// initialize connectvty_mtrx
	connectvty_mtrx.resize(kr*kc);
	for(int i = 0; i < (kr*kc); ++i) {
		connectvty_mtrx[i].resize(kr*kc);
	}
	// Generate 'connectvty_mtrx' based on Mesh
	// topology connection.
	for(int src = 0; src < (kr*kc); ++src) {
		// cout << "Src ["<<src<<"]:\t";
		for(int dst = 0; dst < (kr*kc); ++dst) {
			if (src == dst)
				connectvty_mtrx[src][dst] = -1;
			else if(/*(dst == (src - 1)) ||
					(dst == (src + 1)) ||*/
					(dst == (src + kc)) ||
					(dst == (src - kc))) {
					connectvty_mtrx[src][dst] = 1;
			} else if((src % kc == 0) &&
					(dst == src + 1)) {
					connectvty_mtrx[src][dst] = 1;
			} else if((src % kc == kc - 1) &&
					(dst == src - 1)) {
					connectvty_mtrx[src][dst] = 1;
			} else if( 	(src % kc != 0) &&
						(src % kc != kc - 1) &&
						((dst == (src -1)) ||
						(dst == (src + 1)))) {
					connectvty_mtrx[src][dst] = 1;
			} else
				connectvty_mtrx[src][dst] = 0;
		}
	}
	// initialize the 'visited' matrix...
	visited.resize(kr*kc);
	for(unsigned i = 0; i < visited.size(); ++i) {
		visited[i] = false;
	}
	// initialize the vector of vectors of vectors of vectors for routing table
	all_pairs_all_paths.resize(kr*kc);
	for(int src_= 0; src_ < (kr*kc); ++src_) {
		all_pairs_all_paths[src_].resize(kr*kc);
	}
	return;
};

/*******************************
This function removes link while keeping
 the topology connected
 'new_connectvty_mtrx' will contain the
 new random topology information.
********************************/
// it will do the sanity check...
// and update the new connectivity matrix
bool
gen_topology::remove_link(int src_node) {
	bool success = false;
	int cntr_ = 0;
	vector<int> col_;
	redo:  // TODO: remove goto statement
	col_ = new_connectvty_mtrx[src_node];

	int sum_of_elems = std::accumulate(col_.begin(), col_.end(), 0);

	if (sum_of_elems > 0) {
		unsigned long dest_;
		dest_ = rand() % col_.size();
		while((col_[dest_] == 0) || (col_[dest_] == -1)) {
			// cout << "stuck!!" << endl;
			dest_ = rand() % col_.size();

		}

//		cout << "came-out!!" << endl;
//		cout << "src_node: " << src_node << endl;
//		cout << "dest_: " << dest_ << endl;
		// cout << "col_[dest_]: " << col_[dest_] << endl;
		assert(col_[dest_] == 1);

		// at this point 'dest_' is found connected
		// only remove this if the 'src_' row for this dest has
		// 'sum_of_elems' > 0... otherwise redo: [find new 'dest_']
		std::vector<int> dst_row_;
		dst_row_ = new_connectvty_mtrx[dest_];
		int dst_sum_of_elems = std::accumulate(dst_row_.begin(), dst_row_.end(), 0);

		if (dst_sum_of_elems > 0) {
			assert(new_connectvty_mtrx[src_node][dest_] == 1);
			assert(new_connectvty_mtrx[dest_][src_node] == 1);
			new_connectvty_mtrx[src_node][dest_] = 0;
			new_connectvty_mtrx[dest_][src_node] = 0;
			// you are good to remove these pairs..
			// cout << "-------new_connectvty_mtrx--------" << endl;
			// print_matrix(new_connectvty_mtrx);

			/* Check if graph is 'connected' or not.. if not
			 * then re-do...*/
			success = check_connected();
			if(success == false) {
				/* 1. roll over */
				new_connectvty_mtrx[src_node][dest_] = 1;
				new_connectvty_mtrx[dest_][src_node] = 1;
				/* 2. increment counter */
				cntr_++;
				/* 3. check if limit reached.. return */
				if(cntr_ > 10) {
					return success;
				} else {
					/* 3. goto redo */
					goto redo;
				}
			}
			else if(success == true) {
				/* do nothing */
				/*cout << "Removed link between nodes: " <<
				    src_node << "----" << dest_ << endl;*/
			}
		} else {
			cntr_++;
			if(cntr_ > 10) {
				return success;
			}
			goto redo;
			// try again...
		}
	} else {
		// do nothing
	}

	return success;
}

/*******************************
Describe function
 Initialize the 'new_connectvty_matrx' for a
 new random irregular topology on each call to
 this function.
********************************/
void
gen_topology::generate_topology() {

    // initialize the 'new connectvty mtrx
    // with original one
	new_connectvty_mtrx = connectvty_mtrx;

	// cout << "-------new_connectvty_mtrx--------" << endl;
	// print_matrix(new_connectvty_mtrx);

	for(int link = 0; link < links_removed; ++link) {
		int cntr = 0;
		bool rslt = false;
//		repeat:
		do {
            cntr++;
            cout << "gen_topology::generate_topology():: cntr:" << cntr << endl;
            if (cntr >= 1000) {
                cout << "exhaousted the number of trials; "\
                        "couldn't find link to remove.. exiting\n";
                exit(-1);
            }
            int rand_src_ = rand() % (kr * kc); // generate random num.
            // remove one link of this random node if sum of row >= 0
            // repeat otherwise... with maximum upto 1000 times..
            // otherwise giveup...
            rslt = remove_link(rand_src_);
        }while (rslt == false); // need elegant way to do this...
	}

	return;
}

/*********************************************
 * findAllPathsUtils():
 * *******************************************/
bool
gen_topology::findAllPathsUtil(int src, int dst,
								vector<bool>& visited,
							   vector<int>& path, int &path_index,
							   int &orig_src) {
	// Mark the current node and store it in path
	visited[src] = true;
	path[path_index] = src;
	path_index++;
	// If current vertex is same as destination, then print
	// current_path[]
	if (src == dst) {
        // should I add the path in the routing table?
        bool add = true;
	    vector<int> tmp(&path[0], &path[path_index]);
	    // add only those paths which are up/dn routing compliant.
	    // 'up_down_mtrx'
	    vector<char>valid_path;
	    for(unsigned itr_=0; itr_ < (tmp.size() - 1); ++itr_) {
	        char curr_char = up_down_mtrx[tmp[itr_]][tmp[itr_+1]];
	        if(valid_path.size() > 1) {
	            char prev_char = valid_path.back();
	            if ((prev_char == 'd')
                    && (curr_char == 'u')) {
                    add = false; // don't add it to path
                    assert(0); // this case shouldn't arrieve in first place.
                    // cout << "cannot add because of path violation" << endl;
	            }
	        }
	        valid_path.push_back(curr_char);
	    }
//	    cout << "up/dn link traversal is as follows: " << endl;
//	    for(unsigned j=0; j < valid_path.size(); ++j) {
//	        cout << valid_path[j] << " ";
//	    }
	    // cout << endl;
	    // only add to routing table
	    if(add == true) {
            // all_pairs_all_paths[src_][dst_].push_back(tmp);
            for (unsigned s = 0; s < tmp.size(); ++s) { // 's' and 'd' are indexes
                for (unsigned d = s + 1; d < tmp.size(); ++d) {
                    vector<int> tmp_1(&tmp[s], &tmp[d + 1]);
                    // this means there's already an entry present with given
                    // source-destination pair
                    if (all_pairs_all_paths[tmp_1.front()][tmp_1.back()].size() >= 1)
                        continue;
                    all_pairs_all_paths[tmp_1.front()][tmp_1.back()].push_back(tmp_1);
                }
            }
        }
	    /*
		for (int i = 0; i < path_index; i++) {
			cout << path[i] << " ";
		}
		cout << endl;*/
	}
	else {
		// If current vertex is not destination
		// Recur for all vertices adjacent to current vertex
		for(unsigned itr = 0; itr < new_connectvty_mtrx[src].size(); ++itr) {
			if(new_connectvty_mtrx[src][itr] == 1) {
				if(!visited[itr]) {
                    // you don't allow for recursion if there's 'd'->'u' turn
                    // from the current path--vector
                    if (path_index >= 2) {
                        int prev_node = path[path_index - 2];
                        int cur_node = path[path_index - 1];
                        if ((up_down_mtrx[prev_node][cur_node] == 'd') &&
                            (up_down_mtrx[cur_node][itr] == 'u'))
                            continue; // don't schedule an dn->up turn
                    }
				    if (findAllPathsUtil(itr, dst, visited, path, path_index, orig_src)
				    	== true)
				    	break;
                }
			}
		}
	}
	// Remove current vertex from the path and mark it as unvisited.
	path_index--;
	visited[src] = false;
	// TODO: you can put logic here to recurse back until there's
	// NO 'd' -> 'u' turn
	if(all_pairs_all_paths[orig_src][dst].size() >= 1) {
		return true;
	} else
		return false;
}

/*********************************************
 * findAllPaths():
 * *******************************************/
void
gen_topology::findAllPaths(int src, int dst) {
	// Mark all the vertices as not visited
	vector<bool> visited;
	visited.resize(new_connectvty_mtrx.size());
	// Create an array to store paths
	vector<int> path;
	path.resize(new_connectvty_mtrx.size());
	int path_index = 0; // initialize path[] as empty

	// Initialize all vertices as not visited.
	for (unsigned i = 0; i < new_connectvty_mtrx.size(); i++)
		visited[i] = false;

	// Call the recursive helper function to print all paths
	findAllPathsUtil(src, dst, visited, path, path_index, src);
}

/**************************
 *
 * **************************/
void
gen_topology::populate_all_pair_all_paths() {
    // re-initialize the 'all_pairs_all_paths' to avoid clobbering
    // for differnt connectivity matrix in each run.
    all_pairs_all_paths.clear();
    all_pairs_all_paths.resize(kr*kc);
    for(int src_= 0; src_ < (kr*kc); ++src_) {
        all_pairs_all_paths[src_].resize(kr*kc);
    }
    findAllPaths(3,8); // comment it out
	for(unsigned src = 0; src < new_connectvty_mtrx.size(); ++src) {
		for(unsigned dst = 0; dst < new_connectvty_mtrx.size(); ++dst) {
			if(src != dst) {
				// cout << "All paths between " << src << "---" << dst << endl;
				// Do not schedule this function if there's already entry in
				// the 'all_pairs_all_path' for this source-destination pair
				if(all_pairs_all_paths[src][dst].size() >= 1)
					continue;
				findAllPaths(src, dst);
			}
		}
	}
}

/*******************************
This function initilizes the link with
 'up/dn' based on the irregular topology
 with root node 0.
********************************/
void
gen_topology::populate_up_dwn_mtrx() {
	std::vector< std::vector<int> > shortest_distance;
	shortest_distance = new_connectvty_mtrx;
	// assign infinity weight for 0-value element and 0 to '-1' value element
	for(unsigned i = 0; i < shortest_distance.size(); ++i) {
		for(unsigned j = 0; j < shortest_distance.size(); ++j) {
			if(shortest_distance[i][j] == -1) {
				shortest_distance[i][j] = 0;
			} else if (shortest_distance[i][j] == 0) {
                shortest_distance[i][j] = 1000; // arbitrary high value
			}
		}
	}
    // cout << "-------shortest_distance-matrix----------" << endl;
	// print_matrix(shortest_distance);

	// API for finding shortest distance.
	bool change = true;
	int nodes = new_connectvty_mtrx.size();
	// cout << "nodes: " << nodes << endl;
	while (change) {
		change = false;
		for (int i = 0; i < nodes; i++) { 		// 'src'
			for (int j = 0; j < nodes; j++) { 	// 'dst'
//				std::vector<int> intermediate_;
				int minimum_ = shortest_distance[i][j];
//				int previous_minimum_ = minimum_;
				for (int k = 0; k < nodes; k++) { // 'intermediate node'
					minimum_ = min(minimum_,
								shortest_distance[i][k] + shortest_distance[k][j]);
				}

				if(shortest_distance[i][j] != minimum_) {
					change = true;
					shortest_distance[i][j] = minimum_;
				}
			}
		}
	}
    // cout << "-------shortest_distance-matrix----------" << endl;
    // print_matrix(shortest_distance);
	// assert(0);
	// assign it to the member function
	all_pairs_shortest_path = shortest_distance;
	// Print function
	// print_matrix(shortest_distance);

	// initialize 'up_down_mtrx' here
	up_down_mtrx.resize(kr*kc);
	for(int i = 0; i < (kr*kc); ++i) {
		up_down_mtrx[i].resize(kr*kc);
	}
	for (int i = 0; i < nodes; ++i) {  // Rows (src)
		for (int j = 0; j < nodes; ++j) { // Cols (dest)
			if(i == j || (new_connectvty_mtrx[i][j] == 0)) {
				// both 'src' and 'dest' are same
				up_down_mtrx[i][j] = '-';
			}
			else if(shortest_distance[0][i] >
				shortest_distance[0][j]) {
				up_down_mtrx[i][j] = 'u';
			} else {
				up_down_mtrx[i][j] = 'd';
			}
		}
	}
	// Print function
	// print_char_matrix(up_down_mtrx);

	// assert(0);

	return;
}

/*******************************
Describe function
********************************/
void
get_params(int argc, const char** argv) {

	for (int ii = 1; ii < argc; ii++) {
		if ((string(argv[ii]) == "-row") ||
			(string(argv[ii]) == "-rows")) {
			if (ii < argc - 1) {
	     		  num_row =  atoi(argv[ii+1]);
			  ii += 1;
			}
	    }
	    else if ((string(argv[ii]) == "-col") ||
		    	(string(argv[ii]) == "-cols")) {
	    	if (ii < argc - 1) {
	    		num_col = atoi(argv[ii+1]);
	    		ii += 1;
	    	}
	    }
	    else if ((string(argv[ii]) == "-remove_link") ||
	    		(string(argv[ii]) == "-remove_links")) {
	    	if (ii < argc - 1) {
	    		link_removed = atoi(argv[ii+1]);
	    		ii += 1;
	    	}
	    }
	    else if ((string(argv[ii]) == "-num_topology") ||
	    		(string(argv[ii]) == "-num_topologies")) {
	    	if (ii < argc - 1) {
	    		num_topology = atoi(argv[ii+1]);
	    		ii += 1;
	    	}
	    }
	    else {
			char msg[256];
			printf(msg, "Invalid option %s", argv[ii]);
			usage();
	    }
	}
//	    cout << "num_rows: " << num_row << endl;
//	    cout << "num_cols: " << num_col << endl;
//	    cout << "remove-link: " << link_removed << endl;
//	    cout << "num_topology to be generetated: " << num_topology << endl;
	    // exit(-1);
}

/*******************************
Describe function
********************************/
int main(int argc, char const *argv[])
{

	get_params(argc, argv);
	gen_topology *topology;
	topology = new gen_topology(num_row, num_col, link_removed, num_topology);
	// topology->prnt_state();
	srand((unsigned int) time(NULL));

	// each connectivity matrix is one topology...
	double elapsed_time_ns = 0;
	for(int top_ = 0; top_ < topology->get_num_topology(); ++top_) {
		cout << "====== Call# " << top_ << endl;
		auto start = std::chrono::steady_clock::now(); // start the timer
		topology->generate_topology(); // this is the API generating the topologies...
		auto end = std::chrono::steady_clock::now(); // end the timer
		// find the difference..
		elapsed_time_ns = elapsed_time_ns +
		                    double(std::chrono::duration_cast \
		                            <std::chrono::nanoseconds>(end - start).count());
        int nodes_ = topology->get_kr() * topology->get_kc();
		std::ofstream outFile( std::to_string(nodes_)+"_nodes-connectivity_matrix_"+std::to_string(top_)+
		                        "-links_removed_"+std::to_string(topology->get_links_removed())
		                        +".txt");
		outFile << topology->get_kr() << " " << topology->get_kc() << endl;
		// Do file handling here...
		outFile << "-------Topology--------" << endl;
		// get the start time
		// run some code...
		for(int row = 0; row < (topology->get_kr()*topology->get_kc()); ++row) {
			// outFile << "Row ["<<row<<"]:\t";
			for(int col = 0; col < (topology->get_kr()*topology->get_kc()); ++col) {
				outFile << topology->new_connectvty_mtrx[row][col] << " ";
			}
			outFile << endl;
		}
        outFile << endl;
		// put spinRing input matrix here, derived from topology
		outFile << "-------SpinRing--------" << endl;
		outFile << nodes_ << endl;
		outFile << endl;
        for (unsigned i = 0;
        i < topology->new_connectvty_mtrx.size(); ++i) {
            for (unsigned k = 0;
            k < topology->new_connectvty_mtrx[i].size(); ++k) {
                if(topology->new_connectvty_mtrx[i][k] == 1) {
                    outFile << i << ' ' << k << endl;
                }
            }
        }
		topology->populate_up_dwn_mtrx();
        outFile << endl;

        // Put up/dn link-labels in the same file:
        // Do file handling here...
        outFile << "-------UP/DOWN--------" << endl;
        for(int row = 0;
            row < (topology->get_kr()*topology->get_kc()); ++row) {
            for(int col = 0;
                col < (topology->get_kr()*topology->get_kc()); ++col) {
                outFile << topology->up_down_mtrx[row][col] << " ";
            }
            outFile << endl;
        }
        outFile << endl;
		// We have updated the 'up/dn' matrix for this topology.
		// Now find all-pair all-paths for this topology
		topology->populate_all_pair_all_paths();
		// print out the 'all_pairs_all_paths' structure here
        outFile << "-------UP/DOWN_PATHS--------" << endl;
		for(unsigned src = 0; src < topology->all_pairs_all_paths.size(); ++src) {
			for(unsigned dst = 0; dst < topology->all_pairs_all_paths[src].size(); ++dst) {
				if(src == dst) {
					// do nothing: skip
				}
				else {
					outFile << " [ "<<src<<"]-----:";
					outFile << " ["<< dst << "] :" << endl;
					for(unsigned k=0; k < topology->all_pairs_all_paths[src][dst].size(); ++k) {
						for(unsigned j=0; j < topology->all_pairs_all_paths[src][dst][k].size(); ++j) {
							outFile << topology->all_pairs_all_paths[src][dst][k][j] << " ";
						}
						outFile << endl;
					}
				    outFile << endl;
				}
			}
//			outFile << endl;
		}
		outFile.close();

	}
	// output
	cout << "Elapsed Time (ms): " << elapsed_time_ns / 1e6 << endl;
	return 0;
}


/*******************************
Describe function
********************************/
void
usage() {
	cout << "sample command to run is:" << endl;
	cout << "./gen_topology <rows> <cols> <link-removed> <number-of-topology-gen>" << endl;
	cout << "example: " << endl;
	cout << " ./gen_topology -row 3 -col 5 -remove_links 3 -num_topology 4" << endl;
	cout << "OR" << endl;
	cout << " ./gen_topology -rows 3 -cols 5 -remove_links 3 -num_topologies 4" << endl;
	exit(-1);
	return;
}