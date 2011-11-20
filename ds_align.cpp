#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <dirent.h>
#include <iostream>
#include <fstream>
#include "edit_distance.h" 

using namespace std;


void find_and_replace(string& in_str, const string& original_string, const string& replacement_string)
{
  size_t n = 0;
  while((n = in_str.find(original_string, n)) != std::string::npos) {
     in_str.replace(n, original_string.length(), replacement_string);
     n += replacement_string.length();
  }
}

// Return a vector of strings containing all files call original.tif under data_folder
vector<string> get_all_edgewalks(string data_folder) {

  cout << "Data folder: " << data_folder << endl;
  vector<string> file_list;
  DIR *dp;

  dp = opendir (data_folder.c_str());

  struct dirent *ep;
  while ((ep = readdir (dp)) != NULL) {
    string filename = ep->d_name;
    if((filename != ".") && (filename != "..")) {
      file_list.push_back(data_folder + "/" + filename + "/edge_walks");
    }
  }

  closedir (dp);

  return file_list;
}

class Score {
public:

  Score(int inscore,int inindex,string inname) {
    score = inscore;
    index = inindex;
    name  = inname;
  }

  bool operator<(const Score &other) const {
    if(score < other.score) return true;
    return false;
  }

  bool operator>(const Score &other) const {
    if(score > other.score) return true;
    return false;
  }

  bool operator==(const Score &other) const {
    if(score == other.score) return true;
    return false;
  }

  string str() {
    string strout;

    strout += stringify(score);
    strout += " ";
    strout += name;
    strout += " ";
    return strout;
  }

  int score;
  int index;
  string name;
};


class FragmentInformation {
public:

  FragmentInformation(string filename_in) {
    ifstream input(filename_in.c_str());

    string temp;
    input >> temp;
    input >> leftside;
    input >> temp;
    input >> rightside;
    input >> temp;
    input >> topside;
    input >> temp;
    input >> bottomside;

    filename = filename_in;

 //   best_left_score  = 100000000;
 //   best_right_score = 100000000;

    cout << "filename: " << filename << endl;
    cout << "lsize:" << leftside.size() << endl;
    cout << "rsize:" << rightside.size() << endl;
  }

  int fgcount() {

    int fg =0;
    for(int n=0;n<leftside.size();n++) if(leftside[n] == '1') fg++;
    for(int n=0;n<rightside.size();n++) if(rightside[n] == '1') fg++;

    return fg;
  }

  string filename;

  string leftside;
  string rightside;
  string topside;
  string bottomside;

  void dump_topscores() {
    sort(left_scores.begin(),left_scores.end());
    sort(right_scores.begin(),right_scores.end());
    sort(top_scores.begin(),top_scores.end());
    sort(bottom_scores.begin(),bottom_scores.end());

    cout << "fragment: " << filename << endl;
    cout << "lefttop: ";
    for(int n=left_scores.size()-1;n>(left_scores.size()-20);n--) cout << left_scores[n].str() << " ";
    cout << endl;

    cout << "righttop: ";
    for(int n=right_scores.size()-1;n>(right_scores.size()-20);n--) cout << right_scores[n].str() << " ";
    cout << endl;

    cout << "toptop: ";
    for(int n=top_scores.size()-1;n>(top_scores.size()-20);n--) cout << top_scores[n].str() << " ";
    cout << endl;

    cout << "bottomtop: ";
    for(int n=bottom_scores.size()-1;n>(top_scores.size()-20);n--) cout << bottom_scores[n].str() << " ";
    cout << endl;



    cout << "fragmentend" << endl;
  }


  bool operator<(const FragmentInformation &other) const {
    Score s  = get_top_score();
    Score os = other.get_top_score();

    if(s < os) return true;
    return false;
  }

  bool operator>(const FragmentInformation &other) const {
    Score s  = get_top_score();
    Score os = other.get_top_score();

    if(s > os) return true;
    return false;
  }

  bool operator==(const FragmentInformation &other) const {
    Score s  = get_top_score();
    Score os = other.get_top_score();

    if(s == os) return true;
    return false;
  }

  Score get_top_score() const {

  vector<Score> left_scores_t = left_scores;
  vector<Score> right_scores_t = right_scores;


    sort(left_scores_t.begin() ,left_scores_t.end());
    sort(right_scores_t.begin(),right_scores_t.end());

    if(left_scores_t[left_scores_t.size()-1] > right_scores_t[right_scores_t.size()-1]) return left_scores_t[left_scores_t.size()-1];
                                                                           else return right_scores_t[right_scores_t.size()-1];
  }


  vector<Score> left_scores;
  vector<Score> right_scores;

  vector<Score> top_scores;
  vector<Score> bottom_scores;
};

int find_frag(int id,vector<FragmentInformation> &frags) {

  for(size_t i=0;i<frags.size();i++) {
    size_t n = 0;
    if((n = frags[i].filename.find(stringify(id), n)) != std::string::npos) return i;
  }
  return -1;

}

int main(int argc, char* argv[]) {

  cerr << "ds_align <data folder>" << endl;
  cerr << "ds_align frag1 frag2" << endl;


  cout << " argc: " << argc << endl;
  if(argc == 3) {
    cout << "aligning: " << argv[1] << " " << argv[2] << endl;
    FragmentInformation f1(argv[1]);
    FragmentInformation f2(argv[2]);
   
 /*
    cout << "f1 left trim  : " << trim(f1.leftside)   << endl;
    cout << "f1 right trim : " << trim(f1.rightside)  << endl;
    cout << "f1 top trim   : " << trim(f1.topside)    << endl;
    cout << "f1 bottom trim: " << trim(f1.bottomside) << endl;
    cout << "f2 left trim  : " << trim(f2.leftside)   << endl;
    cout << "f2 right trim : " << trim(f2.rightside)  << endl;
    cout << "f2 top trim   : " << trim(f2.topside)    << endl;
    cout << "f2 bottom trim: " << trim(f2.bottomside) << endl;
*/

    cout << "PROCESS LEFT" << endl;
    int lscore = edit_distance(f1.leftside  ,f2.rightside);
    cout << "PROCESS RIGHT" << endl;
    int rscore = edit_distance(f1.rightside ,f2.leftside);
    cout << "PROCESS TOP" << endl;
    int tscore = edit_distance(f1.topside   ,f2.bottomside);
    cout << "PROCESS BOTTOM" << endl;
    int bscore = edit_distance(f1.bottomside,f2.topside);

    cout << "Left score  : " << lscore << endl;
    cout << "Right score : " << rscore << endl;
    cout << "Top score   : " << tscore << endl;
    cout << "Bottom score: " << bscore << endl;
    return 0;
  }

  // Load data
  vector<string> walk_filenames = get_all_edgewalks(string(argv[1]));

  vector<FragmentInformation> fragments;
  for(int n=0;n<walk_filenames.size();n++) {
    cerr << "Loading: " << walk_filenames[n] << endl;

    string walk_filename = walk_filenames[n];

    fragments.push_back(FragmentInformation(walk_filename));
  }








/*
  size_t only_fragment = find_frag(142,fragments);
  bool only_fragment_flag = true;
  if(only_fragment_flag == true) {
    // Find best scores for each fragment (currently left and right only)
    #pragma omp parallel for
    for(int j=0;j<fragments.size();j++) {
      int i = only_fragment;
      cout << "aligning: " << i << " " << j << endl;
      int lscore = edit_distance(fragments[i].leftside,fragments[j].rightside);
      int rscore = edit_distance(fragments[i].rightside,fragments[j].leftside);
      //int tscore = edit_distance(fragments[i].topside,fragments[j].bottomside);
      //int bscore = edit_distance(fragments[i].bottomside,fragments[j].topside);


      #pragma omp critical 
      {
        fragments[i].left_scores  .push_back(Score(lscore,j,fragments[j].filename));
        fragments[i].right_scores .push_back(Score(rscore,j,fragments[j].filename));
        //fragments[i].top_scores   .push_back(Score(tscore,j,fragments[j].filename));
        //fragments[i].bottom_scores.push_back(Score(bscore,j,fragments[j].filename));
      }
    }

    // Dump fragments
    fragments[only_fragment].dump_topscores();
    return 0;
  }
*/


  // Process normally.




  // Find best scores for each fragment (currently left and right only)
  #pragma omp parallel for
  for(int i=0;i<fragments.size();i++) {
    for(int j=0;j<fragments.size();j++) {
      cout << "aligning: " << i << " " << j << endl;
      int lscore = edit_distance(fragments[i].leftside,fragments[j].rightside);
      fragments[i].left_scores.push_back(Score(lscore,j,fragments[j].filename));

      int rscore = edit_distance(fragments[i].rightside,fragments[j].leftside);
      fragments[i].right_scores.push_back(Score(rscore,j,fragments[j].filename));

      int tscore = edit_distance(fragments[i].topside,fragments[j].bottomside);
      fragments[i].top_scores.push_back(Score(tscore,j,fragments[j].filename));
 
      int bscore = edit_distance(fragments[i].bottomside,fragments[j].topside);
      fragments[i].bottom_scores.push_back(Score(bscore,j,fragments[j].filename));
    }
  }

  sort(fragments.begin(),fragments.end());

  // Dump fragments
  for(size_t n=0;n<fragments.size();n++) {
    fragments[n].dump_topscores();
  }

}
