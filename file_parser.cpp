//
// Created by Mayank Parasar on 2018-11-09.
//
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;
// Globals
int xlen, ylen;

struct entry {
    int next_router_id;
    string direction_;
    entry() :
        next_router_id(-1),
        direction_("Unknown")
    {}
    entry(int id, string dirn) :
            next_router_id(id),
            direction_(dirn)
    {}
};

struct upDn_ {
    int src;
    int dst;

    bool operator==(const upDn_ &pair_) const {
        return (src == pair_.src && dst == pair_.dst);
    }

    bool operator<(const upDn_ &pair_)  const {
        return ((src < pair_.src) ||
                (src == pair_.src && dst < pair_.dst));
    }

    upDn_(int src_, int dst_) :
        src(src_), dst(dst_)
    {}
};
vector<vector< vector< entry> > > routingTable;

void
populate_routingTable(std::vector<int>& path_) {
    int dst = path_.back();
    int src = path_.front();
    entry entry_;
    for(int curr_ = 0, nxt_ = curr_ + 1;
        curr_ < path_.size() && nxt_ < path_.size();
        curr_++, nxt_++) {
        string dirn_;
        if(path_[nxt_] == (path_[curr_] - 1)) {
            // West
            entry_ = {path_[nxt_], "West"};
        }
        else if(path_[nxt_] == (path_[curr_] + 1)) {
            // East
            entry_ = {path_[nxt_], "East"};
        }
        else if(path_[nxt_] == (path_[curr_] + ylen)) {
            // North
            entry_ = {path_[nxt_], "North"};
        }
        else if(path_[nxt_] == (path_[curr_] - ylen)) {
            // South
            entry_ = {path_[nxt_], "South"};
        }
        else if(path_[nxt_] == path_[curr_]) {
            // skip do nothing...
        }
        else {
            cout << " this is not possible" << endl;
            assert(0);
        }

        // push the entry_ into routingTable
        // only push if entry is unique
        // if(routingTable[path_[curr_]][dst].size() == 0)
            // routingTable[path_[curr_]][dst].push_back(entry_);
        // New logic is to put the complete path for given
        // src and destination pair
        routingTable[src][dst].push_back(entry_);
    }
}

int main(int argc, char const *argv[]) {
    ifstream inFile(argv[1]);
    string word;

    inFile >> word;
    xlen = stoi(word);  // rows
    inFile >> word;
    ylen = stoi(word);  // cols

    // Resize the table
    routingTable.resize(xlen*ylen);
    for(int i = 0; i < xlen*ylen; ++i) {
        routingTable[i].resize(xlen*ylen);
    }
    bool top_ = false;
    bool spinRing = false;
    bool up_dn = false;
    bool up_dn_path = false;
    bool path_start = false;
    bool path_end = false;
    // bool false_positive = true;
    std::vector<int> tmp_path;

    while(!(inFile.eof())) {
        inFile >> word;

        if((word.find("Topology") != -1)) {
            top_ = true;
            spinRing = false;
            up_dn = false;
            up_dn_path = false;
            // cout << "got the word topology" << endl;
        }
        if((word.find("SpinRing") != -1)) {
            top_ = false;
            spinRing = true;
            up_dn = false;
            up_dn_path = false;
        }
        if((word.find("UP/DOWN") != -1)) {
            top_ = false;
            spinRing = false;
            up_dn = true;
            up_dn_path = false;
        }
        if((word.find("UP/DOWN_PATHS") != -1)) {
            top_ = false;
            spinRing = false;
            up_dn = false;
            up_dn_path = true;
        }

        /*Based on flag value do the work*/
        if( top_ == true ) {
            // cout << word << endl;
            /*Do work on this subsection*/
        }
        /*********************************/
        /*if( up_dn == true &&
            up_dn_path == false) {
            // populate the map here
            cout << word << endl;
            cout.flush();
            if(word == "d" ||
                word == "u") {
                cout << word;
                // cout.flush();
            }


        } else*/
        if( up_dn_path == true ) {

            if (inFile.peek() == EOF) {
                path_start = false;
                path_end = true;
            }
            if((path_start == false) &&
               (path_end == true) &&
                    (tmp_path.size()>0)) {
                if(tmp_path[0] == 63 && tmp_path.back() == 54)
                    cout << "Interesting case" << endl;
                populate_routingTable(tmp_path);
                for (int idx = 0; idx < tmp_path.size(); ++idx) {
                    cout << tmp_path[idx] << " ";
                }
                cout << endl;
                tmp_path.clear();
                cout << "tmp_path size after clear(): "\
                        << tmp_path.size();
                cout << endl;

            }
            if (word =="[" || inFile.peek() == EOF) {
                path_start = false;
                path_end = true;
            }
            if (path_start == true &&
                path_end == false) {
                // cout << stoi(word);
                tmp_path.push_back(stoi(word));
                if(tmp_path[0] == 63 && tmp_path.size() >= 26)
                    cout << "Interesting case" << endl;
            }
            if (word == ":") {
                path_start = true;
                path_end = false;
            }

            // cout << "reached endof of the line" << endl;
            // cout << word << endl;
        }
    }

    // re-read the file..
    // this time line by line..
    inFile.clear(); // this was important
    inFile.seekg(0, std::ios::beg);
    string line;
    int src = 0;
    int dst = 0;
    bool start = false;
    bool stop = false;
    map<upDn_, char> global_upDn;
    while (std::getline(inFile, line)) {

        if( start == true &&
            line.empty()) {
            start = false;
            stop = true;
        }

        if( start == true &&
            stop == false) {
            cout << line << endl;
            //break this line into deliminter
            for(auto x : line) {

                if(x == 'u') {
                    cout << x << endl;
                    pair<upDn_, char> p((upDn_{src,dst}),x);
                    global_upDn.insert(p);
                }
                if(x == 'd') {
                    cout << x << endl;
                    pair<upDn_, char> p((upDn_{src,dst}),x);
                    global_upDn.insert(p);
                }
                if(x == ' ') {
                    // do not increment dst here here
                } else {
                    dst += 1;
                }
            }
            dst = 0; // reset
            src += 1; // increment.
            cout.flush();
        }

        if((line.find("UP/DOWN") != -1)) {
            // cout << line << endl;
            cout. flush();
            start = true;
        }



    }

    inFile.close();

    // cout global--map
    for (auto& t : global_upDn)
        std::cout << t.first.src << " "
                  << t.first.dst << " "
                  << t.second << " "
                  << "\n";
    // check a certain value:
    upDn_ tmp(5,4);
//    cout << "tmp(5,4): " << global_upDn[tmp] << endl;
    cout << "tmp(5,4): " << global_upDn.at(tmp) << endl;
    tmp = {5,5};
    // cout << "tmp(5,5): " << global_upDn.at(tmp) << endl;
    // cout routing table:
    for(int row = xlen*ylen - 1; row >= 0; --row ) {
        cout << "Row ["<< row <<"]:\t";
        for(int col = 0; col < xlen*ylen; ++col) {
            cout << "(";
            for (auto i : routingTable[row][col]) {
                cout << i.direction_ <<", " << i.next_router_id;
                cout << " ; ";
            }
            cout << ") ";
        }
        cout << endl;
    }
	return 0;
}