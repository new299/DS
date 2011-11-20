#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {

  cerr << "ds_merge image1.tif image2.tif outputfile.tif ox oy" << endl;
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
    //i1.remove_fragments(0x88888888,0,30);
    //i1.blockify();
    //i1.remove_fragments(0x88888888,0,110);
    i1 = i1.erode(0);
    i1 = i1.erode(0);
    i1.linear_repair(0,0x88888888,2);
    i1.remove_fragments(0x88888888,0,2);
    i1.linear_repair(0,0x88888888,6);
    i1.remove_fragments(0x88888888,0,10);
    i1.linear_repair(0,0x88888888,10);
    i1.remove_fragments(0x88888888,0,30);

    i1 = i1.add_border(80,0);
 
    Image i2;
    i2.load_tiff(input2_filename);
    //i2.load_tiff_8bit_grey(input2_filename);
    i2.make_greyscale();
    i2.threshold(0.50,0,0x88888888,0xFFFFFFFF);
    //i2.remove_fragments(0x88888888,0,30);
    //i2.blockify();
    //i2.remove_fragments(0x88888888,0,110);
    i2 = i2.erode(0);
    i2 = i2.erode(0);
    i2.linear_repair(0,0x88888888,2);
    i2.remove_fragments(0x88888888,0,2);
    i2.linear_repair(0,0x88888888,6);
    i2.remove_fragments(0x88888888,0,10);
    i2.linear_repair(0,0x88888888,10);
    i2.remove_fragments(0x88888888,0,30);

    int quality;
    int ox = convertTo<int>(argv[4]);
    int oy = convertTo<int>(argv[5]);

    quality = i1.get_merge_score(i2,ox,oy,0);
    cout << "Merge quality: " << quality << endl;
    cout << "Merge position: " << ox << " " << oy << endl;
 
    i1.overlay(i2,ox,oy,0);
    
    i1.save_tiff_grey_8bit(output_filename);

    return 0;
  }

}
