#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <dirent.h>

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
vector<string> get_all_original(string data_folder) {

  vector<string> file_list;
  DIR *dp;

  dp = opendir (data_folder.c_str());

  struct dirent *ep;
  while ((ep = readdir (dp)) != NULL) {
    string filename = ep->d_name;
    if((filename != ".") && (filename != "..")) {
      file_list.push_back(data_folder + "/" + filename + "/original.tif");
    }
  }

  closedir (dp);

  return file_list;
}

int main(int argc, char* argv[]) {

  cerr << "ds_orientate <data folder> [outputfile]" << endl;
  cerr << "If you supply a single tif file instead of <data folder> you must also supply outputfile" << endl;

  string arg = argv[1];
  if(arg.substr(arg.length()-3,3) == "tif") {
    cerr << "Processing single file" << endl;
    Image i;
    i.load_tiff(arg);
    Image i2 = i.orientate(Image::frag_type_arrow,0);
    i2.save_tiff_rgb(string(argv[2]));
    return 0;
  }

  vector<string> image_filenames = get_all_original(string(argv[1]));

  for(size_t n=0;n<image_filenames.size();n++) {
    cerr << "Processing: " << image_filenames[n] << endl;
    Image i;
    i.load_tiff(image_filenames[n]);
    Image i2 = i.orientate(Image::frag_type_arrow,0);
    string output_filename = image_filenames[n];
    find_and_replace(output_filename,"original","orientated");
    i2.save_tiff_rgb(output_filename);
  }

}
