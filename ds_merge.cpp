#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {

  cerr << "ds_merge image1.tif image2.tif outputfile.tif" << endl;
  cerr << "If you supply a single tif file instead of <data folder> you must also supply outputfile" << endl;

  string arg = argv[1];
  if(arg.substr(arg.length()-3,3) == "tif") {

    string input1_filename = argv[1];
    string input2_filename = argv[2];
    string output_filename = argv[3];

    Image i1;
    //i1.load_tiff_8bit_grey(input1_filename);
    i1.load_tiff(input1_filename);
    i1.make_greyscale();
    i1.threshold(0.5,0,0x88888888,0xFFFFFFFF);
    i1.linear_repair(0,0x88888888);
    i1.remove_fragments(0x88888888,0,30);
    i1.blockify();
    i1.remove_fragments(0x88888888,0,110);

    Image i2;
    i2.load_tiff(input2_filename);
    //i2.load_tiff_8bit_grey(input2_filename);
    i2.make_greyscale();
    i2.threshold(0.50,0,0x88888888,0xFFFFFFFF);
    i2.linear_repair(0,0x88888888);
    i2.remove_fragments(0x88888888,0,30);
    i2.blockify();
    i2.remove_fragments(0x88888888,0,110);


    int quality;
    int ox;
    int oy;
    Image nimage = i1.merge(i2,0,quality,ox,oy);
    cout << "Merge quality: " << quality << endl;
    cout << "Merge position: " << ox << " " << oy << endl;
 
    cout << "nimage size: " << nimage.m_width << " " << nimage.m_height << endl;
    nimage.save_tiff_grey_8bit(output_filename);

    return 0;
  }

}
