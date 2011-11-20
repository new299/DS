#include "Image.h"
using namespace std;

int main(int argc, char* argv[]) {

  cout << "ds_isolate <input image> <output directory> <subtracted image>" << endl;
cout << "argc: " << argc << endl;
  string input_image_filename = argv[1];
  string output_directory     = argv[2];
  string subtracted_image_filename;
  if(argc > 3) subtracted_image_filename = argv[3];


  Image input;

  input.load_tiff(input_image_filename);

  map<uint32_t,int> pixels1 = input.random_walk_collect(20              ,               20,50000000,10);
  map<uint32_t,int> pixels2 = input.random_walk_collect(input.m_width-20,               20,50000000,10);
  map<uint32_t,int> pixels3 = input.random_walk_collect(              20,input.m_height-20,50000000,10);
  map<uint32_t,int> pixels4 = input.random_walk_collect(input.m_width-20,input.m_height-20,50000000,10);

  map<uint32_t,int> pixels; // = input.random_walk_collect(20,20,1000000000,10);
  pixels.insert(pixels1.begin(),pixels1.end());
  pixels.insert(pixels2.begin(),pixels2.end());
  pixels.insert(pixels3.begin(),pixels3.end());
  pixels.insert(pixels4.begin(),pixels4.end());

  cerr << "Acquired " << pixels.size() << " during random walk" << endl;

  pixels = grow_pixel_set(pixels,10);
  cerr << "Grew pixel set to: " << pixels.size() << endl;
  input.set_pixels(pixels,0,0);
  input.remove_isolated_pixels(0);

  if(subtracted_image_filename.length() != 0) input.save_tiff_rgb(subtracted_image_filename);
  cerr << "Cleared background pixels, performing segmentation..." << endl;

  vector<Image> fragments = input.segment_on_value(0);

  // Create output directory folder
  mkdir(output_directory.c_str(),0777);

  cerr << "Saving fragments..." << endl;
  for(size_t n=0;n<fragments.size();n++) {
    if(fragments[n].pixel_count() > 20) {
      string frag_dir = output_directory + string("/") + string("fragment_") + stringify(n);

      mkdir(frag_dir.c_str(),0777);
      fragments[n].save_tiff_rgb(frag_dir + string("/") + string("original.tif"));
    }
  }

}
