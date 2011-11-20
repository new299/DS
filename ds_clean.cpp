#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <dirent.h>
#include <iostream>
#include <fstream>

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
vector<string> get_all_orientated(string data_folder) {

  vector<string> file_list;
  DIR *dp;

  dp = opendir (data_folder.c_str());

  struct dirent *ep;
  while ((ep = readdir (dp)) != NULL) {
    string filename = ep->d_name;
    if((filename != ".") && (filename != "..")) {
      file_list.push_back(data_folder + "/" + filename + "/orientated.tif");
    }
  }

  closedir (dp);

  return file_list;
}

void process_clean(string input_filename,string output_filename,string walk_filename) {
  cerr << "Input file : " <<  input_filename << endl;
  cerr << "Output file: " <<  output_filename << endl;
  cerr << "Walk  file : " <<  walk_filename << endl;
  Image i;
  i.load_tiff(input_filename);
  Image i2 = i.clean(0);
  i2.make_greyscale();
  i2.fill_hole(0);
  i2.threshold(0.57,0,0x88888888,0xFFFFFFFF);
  i2.remove_isolated_pixels(0);
  i2.linear_repair(0xFFFFFFFF,0x88888888);
  i2.remove_fragments(0x88888888,0);
  i2.remove_hangs(0);
  i2.remove_isolated_pixels(0);
  i2.remove_hangs(0);
  i2.save_tiff_grey_8bit(output_filename);


  size_t topleft_x    ,topleft_y;
  size_t bottomleft_x ,bottomleft_y;
  size_t topright_x   ,topright_y;
  size_t bottomright_x,bottomright_y;

  i2.locate_top_left    (topleft_x    ,topleft_y,0);
  i2.locate_bottom_left (bottomleft_x ,bottomleft_y,0);
  i2.locate_top_right   (topright_x   ,topright_y,0);
  i2.locate_bottom_right(bottomright_x,bottomright_y,0);

  cout << "top left : " << topleft_x     << " " << topleft_y << endl;
  cout << "top right: " << topright_x    << " " << topright_y << endl;
  cout << "bot left : " << bottomleft_x  << " " << bottomleft_y << endl;
  cout << "bot right: " << bottomright_x << " " << bottomright_y << endl;

  ofstream outputwalks(walk_filename.c_str());

  // walk left side
  cout << "walking leftside" << endl;
  outputwalks << "LW ";
  vector<uint32_t> leftside_walk_values;
  vector<Position> leftside_walk_positions;
  bool ret1 = i2.walk_pixels(topleft_x,topleft_y,bottomleft_x,bottomleft_y,0,0x88888888,false,leftside_walk_values,leftside_walk_positions);
  if(ret1) for(size_t n=0;n<leftside_walk_values.size();n++) {if(leftside_walk_values[n] == 0x88888888) outputwalks << "1"; else outputwalks << "0";}
      else outputwalks << "0";

  outputwalks << endl;

  // walk right side - clockwise
  cout << "walking rightside" << endl;
  outputwalks << "RW ";
  vector<uint32_t> rightside_walk_values;
  vector<Position> rightside_walk_positions;
  bool ret2 = i2.walk_pixels(topright_x,topright_y,bottomright_x,bottomright_y,0,0x88888888,true,rightside_walk_values,rightside_walk_positions);
  if(ret2) for(size_t n=0;n<rightside_walk_values.size();n++) {if(rightside_walk_values[n] == 0x88888888) outputwalks << "1"; else outputwalks << "0";}
      else outputwalks << "0";
  outputwalks << endl;

  // walk top - clockwise
  cout << "walking topside" << endl;
  outputwalks << "TW ";
  vector<uint32_t> topside_walk_values;
  vector<Position> topside_walk_positions;
  bool ret3 = i2.walk_pixels(topleft_x,topleft_y,topright_x,topright_y,0,0x88888888,true,topside_walk_values,topside_walk_positions);
  if(ret3) for(size_t n=0;n<topside_walk_values.size();n++) {if(topside_walk_values[n] == 0x88888888) outputwalks << "1"; else outputwalks << "0";}
      else outputwalks << "0";
  outputwalks << endl;

  // walk bottom - anticlock
  cout << "walking bottomside" << endl;
  outputwalks << "BW ";
  vector<uint32_t> bottomside_walk_values;
  vector<Position> bottomside_walk_positions;
  bool ret4 = i2.walk_pixels(bottomleft_x,bottomleft_y,bottomright_x,bottomright_y,0,0x88888888,false,bottomside_walk_values,bottomside_walk_positions);
  if(ret4) for(size_t n=0;n<bottomside_walk_values.size();n++) {if(bottomside_walk_values[n] == 0x88888888) outputwalks << "1"; else outputwalks << "0";}
      else outputwalks << "0";
  outputwalks << endl;

  if((ret1 == false) || 
     (ret2 == false) ||
     (ret3 == false) ||
     (ret4 == false)) {

    cerr << "****************************************** FAILED" << input_filename << endl;
  }
}

int main(int argc, char* argv[]) {

  cerr << "ds_clean <data folder> [outputfile] [walkfilename]" << endl;
  cerr << "If you supply a single tif file instead of <data folder> you must also supply outputfile" << endl;

  string arg = argv[1];
  if(arg.substr(arg.length()-3,3) == "tif") {

    string input_filename  = argv[1];
    string output_filename = argv[2];
    string walk_filename   = argv[3];

    process_clean(input_filename,output_filename,walk_filename);

    return 0;
  }

  vector<string> image_filenames = get_all_orientated(string(argv[1]));

  #pragma omp parallel for
  for(int n=0;n<image_filenames.size();n++) {
    cerr << "Processing: " << image_filenames[n] << endl;

    string output_filename = image_filenames[n];
    find_and_replace(output_filename,"orientated","clean");

    string walk_filename = image_filenames[n];
    find_and_replace(walk_filename,"orientated.tif","edge_walks");

    process_clean(image_filenames[n],output_filename,walk_filename);
  }

}
