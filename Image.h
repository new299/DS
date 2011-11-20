#ifndef DS_IMAGE_H
#define DS_IMAGE_H

#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include "stringify.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>
#include <algorithm>

#define TIFFSetR(pixel, x) ((unsigned char *)pixel)[0] = x
#define TIFFSetG(pixel, x) ((unsigned char *)pixel)[1] = x
#define TIFFSetB(pixel, x) ((unsigned char *)pixel)[2] = x
#define TIFFSetA(pixel, x) ((unsigned char *)pixel)[3] = x

using namespace std;

double distance(double x1,double y1,double z1,
                double x2,double y2,double z2) {

  double x_d = (x1-x2)*(x1-x2);
  double y_d = (y1-y2)*(y1-y2);
  double z_d = (z1-z2)*(z1-z2);
  double distance = sqrt(x_d+y_d+z_d);
  return distance;
}

class Position {
  public:

    Position() {}

    Position(size_t x,size_t y): m_x(x), m_y(y) {}

    bool operator==(Position &other) {
      if ((m_x == other.m_x) && (m_y == other.m_y)) return true;
      return false;
    }

   bool operator<(const Position &other) const {
     if(m_x < other.m_x) return true;
     if(m_x > other.m_x) return false;

     if(m_y < other.m_y) return true;
     if(m_y > other.m_y) return false;

     return false;
   }

    bool operator!=(Position &other) {
      return !(*this == other);
    }

    string get_str() {
      return stringify(m_x) + " " + stringify(m_y);
    }
     
    size_t m_x;
    size_t m_y;
};

void rotate_point(double x,
                  double y,
                  double o_x,
                  double o_y,
                  double radians,
                  double &r_x,
                  double &r_y) {

  //1. calculate theta
  //tan theta = opp (y) / adj (x)

  double dx = x - o_x;
  double dy = y - o_y;

  if(dx < 0) dx = 0-dx;
  if(dy < 0) dy = 0-dy;

  int q=1;
  if((x >  o_x) && (y >= o_y)) q = 1;
  if((x >= o_x) && (y <  o_y)) q = 2;
  if((x <  o_x) && (y <= o_y)) q = 3;
  if((x <= o_x) && (y >  o_y)) q = 4;
  double theta;
  if(q == 1) theta  = atan(dy/dx);
  if(q == 2) theta  = atan(dx/dy);
  if(q == 3) theta  = atan(dy/dx);
  if(q == 4) theta  = atan(dx/dy);

  double     h = sqrt(dx*dx + dy*dy);

  if((dy==0) || (dx==0)) theta = 0;

  theta -= radians;

  for(;theta < 0;) { theta = 0 - theta; theta = ((2*3.14)/4) - theta; q++; if(q==5) q=1;}

  for(;theta > ((2*3.14)/4);) {theta -= ((2*3.14)/4); q--; if(q==0) q=4; }


  if(q == 1) { r_y = (sin(theta) * h); r_x = (cos(theta) * h); }
  if(q == 2) { r_x = (sin(theta) * h); r_y = (cos(theta) * h); }
  if(q == 3) { r_y = (sin(theta) * h); r_x = (cos(theta) * h); }
  if(q == 4) { r_x = (sin(theta) * h); r_y = (cos(theta) * h); }

  if(q == 1) {r_x = o_x + r_x; r_y = o_y + r_y;}
  if(q == 2) {r_x = o_x + r_x; r_y = o_y - r_y;}
  if(q == 3) {r_x = o_x - r_x; r_y = o_y - r_y;}
  if(q == 4) {r_x = o_x - r_x; r_y = o_y + r_y;}

}

map<uint32_t,int> grow_pixel_set(map<uint32_t,int> pixels,int delta) {

  map<uint32_t,int> new_pixels;

  for(map<uint32_t,int>::iterator i=pixels.begin();i != pixels.end();i++) {

    int r1 = TIFFGetR(i->first);
    int g1 = TIFFGetG(i->first);
    int b1 = TIFFGetB(i->first);

    for(int r_d=0-delta;r_d<(delta+1);r_d++) {
      for(int g_d=0-delta;g_d<(delta+1);g_d++) {
        for(int b_d=0-delta;b_d<(delta+1);b_d++) {

          if(((r1+r_d) >= 0) && ((r1+r_d) < 256)  &&
             ((g1+g_d) >= 0) && ((g1+g_d) < 256)  &&
             ((b1+b_d) >= 0) && ((b1+b_d) < 256)) {
            uint32_t v = 0xFFFFFFFF;
            TIFFSetR(&v,(uint8_t) r1+r_d);
            TIFFSetG(&v,(uint8_t) g1+g_d);
            TIFFSetB(&v,(uint8_t) b1+b_d);

            new_pixels[v]++;
          }
        }
      }
    }
  }

  return new_pixels;
}

class Image {

public:

  Image() {
  }

  enum fragment_type { frag_type_arrow, frag_type_rectangle } ;

  void load_tiff(string input_filename) {

    m_image_data.clear();

    TIFF* tif = TIFFOpen(input_filename.c_str(), "r");
    if (tif) {
      uint32 w, h;
      size_t npixels;
      uint32* raster;

      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
      npixels = w * h;
      m_width = w;
      m_height = h;

      raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
      if (raster != NULL) {
	if (TIFFReadRGBAImageOriented(tif, w, h, raster,ORIENTATION_TOPLEFT, 0)) {
	  for(size_t n=0;n<npixels;n++) m_image_data.push_back(raster[n]);
	}
	_TIFFfree(raster);
      }
      TIFFClose(tif);
    }
  }

  void load_tiff_8bit_grey(string input_filename) {
    m_image_data.clear();

    TIFF* tif = TIFFOpen(input_filename.c_str(), "r");
    if (tif) {
      uint32 w, h;
      size_t npixels;
      uint32* raster;

      TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
      TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
      npixels = w * h;
      m_width = w;
      m_height = h;

      raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
      if (raster != NULL) {
	if (TIFFReadRGBAImageOriented(tif, w, h, raster,ORIENTATION_TOPLEFT, 0)) {
	  for(size_t n=0;n<npixels;n++) {
            //m_image_data.push_back((255-raster[n])*(256*256*256));
            m_image_data.push_back((raster[n])*(256*256*256));
          }
	}
	_TIFFfree(raster);
      }
      TIFFClose(tif);
    }
  }

  void save_tiff_rgb(string output_filename) {
    TIFF *output_image;

    // Open the TIFF file
    if((output_image = TIFFOpen(output_filename.c_str(), "w")) == NULL){
      cerr << "Unable to write tif file: " << output_filename << endl;
    }

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, m_width);
    TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, m_height);
    TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    // Write the information to the file
    TIFFWriteEncodedStrip(output_image, 0, &m_image_data[0], m_width*m_height * 4);

    // Close the file
    TIFFClose(output_image);
  }

  void save_tiff_grey_32bit(string output_filename) {
    TIFF *output_image;

    // Open the TIFF file
    if((output_image = TIFFOpen(output_filename.c_str(), "w")) == NULL){
      cerr << "Unable to write tif file: " << output_filename << endl;
    }

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, m_width);
    TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, m_height);
    TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);

    // Write the information to the file
    TIFFWriteEncodedStrip(output_image, 0, &m_image_data[0], m_width*m_height * 4);

    // Close the file
    TIFFClose(output_image);
  }

  void save_tiff_grey_8bit(string output_filename) {
    TIFF *output_image;

    // Open the TIFF file
    if((output_image = TIFFOpen(output_filename.c_str(), "w")) == NULL){
      cerr << "Unable to write tif file: " << output_filename << endl;
    }

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, m_width);
    TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, m_height);
    TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISWHITE);

    // convert data to 8bit
    vector<uint8_t> data;
    for(size_t n=0;n<m_image_data.size();n++) {
      data.push_back(255-(m_image_data[n]/(256*256*256)));
    }

    // Write the information to the file
    TIFFWriteEncodedStrip(output_image, 0, &data[0], m_width*m_height);

    // Close the file
    TIFFClose(output_image);
  }

  void make_greyscale(uint32_t bg_val=0) {

    for(size_t n=0;n<m_image_data.size();n++) {

      double r = TIFFGetR(m_image_data[n]);
      double g = TIFFGetG(m_image_data[n]);
      double b = TIFFGetB(m_image_data[n]);

      if(m_image_data[n] != bg_val) {
        double grey = (0.3*r) + (0.59*g) + (0.11*b); // See http://en.wikipedia.org/wiki/Grayscale
        if(grey == bg_val) grey++;
        m_image_data[n] = grey * 256 * 256 * 256;
      }
    }
  }

  bool is_on_image(size_t x,size_t y) {
    if(x < 0) return false;
    if(y < 0) return false;
    if((x < m_width) && (y < m_height)) return true;
    return false;
  }

  uint32_t get(size_t x,size_t y) {
    return m_image_data[(y*m_width)+x];
  }

  uint32_t get(Position p) {
    return get(p.m_x,p.m_y);
  }

  void set(size_t x,size_t y,uint32_t value) {
    m_image_data[(y*m_width)+x] = value;
  }

  void set(Position p,uint32_t value) {
    set(p.m_x,p.m_y,value);
  }

  map<uint32_t,int> random_walk_collect(size_t x, ///< Starting x position
                                        size_t y, ///< Starting y position
                                        size_t walk_distance, ///< How long should we continue to walk.
                                        double max_delta ///< Maximum delta while performing walk.
                                       ) {

    // This method performs a random walk around the image from a given starting point.
    // It collects pixels as it goes, so as to extract a collection of background values.

    map<uint32_t,int> pixels_found;

    for(size_t n=0;n<walk_distance;n++) {

      int x_delta;
      int y_delta;

      x_delta = ((int)(rand()%3)) - 1;
      y_delta = ((int)(rand()%3)) - 1;

      for(;!is_on_image(x+x_delta,y+y_delta);) {
        x_delta = ((int)(rand()%3)) - 1;
        y_delta = ((int)(rand()%3)) - 1;
      }


      double r1 = TIFFGetR(get(x,y));
      double g1 = TIFFGetG(get(x,y));
      double b1 = TIFFGetB(get(x,y));

      double r2 = TIFFGetR(get(x+x_delta,y+y_delta));
      double g2 = TIFFGetG(get(x+x_delta,y+y_delta));
      double b2 = TIFFGetB(get(x+x_delta,y+y_delta));

      double d = distance(r1,g1,b1,r2,g2,b2);

      if(d < max_delta) {
        pixels_found[get(x+x_delta,y+y_delta)]++;
        x = x+x_delta;
        y = y+y_delta;
      }
    }

    return pixels_found;
  }

  void set_pixels(map<uint32_t,int> values,uint32_t set_val,double delta) {
    for(size_t n=0;n<m_image_data.size();n++) {
      if(values.find(m_image_data[n]) != values.end()) m_image_data[n] = set_val; else {
      
        if(delta != 0)
        for(map<uint32_t,int>::iterator i=values.begin();i != values.end();i++) {
          double d = distance(TIFFGetR(i->first),TIFFGetG(i->first),TIFFGetB(i->first),
                              TIFFGetR(m_image_data[n]),TIFFGetG(m_image_data[n]),TIFFGetB(m_image_data[n]));
          if(d < delta) m_image_data[n] = set_val;
        }
      }
    }
  }

  void not_set_pixels(map<uint32_t,int> values,uint32_t set_val,double delta) {
    for(size_t n=0;n<m_image_data.size();n++) {
      if(values.find(m_image_data[n]) == values.end()) m_image_data[n] = set_val; else {
      
        if(delta != 0)
        for(map<uint32_t,int>::iterator i=values.begin();i != values.end();i++) {
          double d = distance(TIFFGetR(i->first),TIFFGetG(i->first),TIFFGetB(i->first),
                              TIFFGetR(m_image_data[n]),TIFFGetG(m_image_data[n]),TIFFGetB(m_image_data[n]));
          if(d > delta) m_image_data[n] = set_val;
        }
      }
    }
  }



  vector<Position> get_all_adjacent(size_t x,size_t y,uint32_t bg_val) {

    vector<Position> adjacent_positions;
    for(int d_x=-1;d_x<2;d_x++) {
      for(int d_y=-1;d_y<2;d_y++) {
      
        if(!((d_x == 0) && (d_y == 0))) {

          if(is_on_image(x+d_x,y+d_y)) {
            if(get(x+d_x,y+d_y) != bg_val) { adjacent_positions.push_back(Position(x+d_x,y+d_y)); }
          }
        }
      }
    }

    return adjacent_positions;
  }


  // gets all adjacent pixels but NOT those diagonally connected
  vector<Position> get_all_adjacent_nodiag_val(size_t x,size_t y,uint32_t bg_val,uint32_t frag_val) {

    vector<Position> adjacent_positions;

    if(is_on_image(x+1,y  ) && (get(x+1,y  ) == frag_val)) adjacent_positions.push_back(Position(x+1,y  ));
    if(is_on_image(x-1,y  ) && (get(x-1,y  ) == frag_val)) adjacent_positions.push_back(Position(x-1,y  ));
    if(is_on_image(x  ,y+1) && (get(x  ,y+1) == frag_val)) adjacent_positions.push_back(Position(x  ,y+1));
    if(is_on_image(x  ,y-1) && (get(x  ,y-1) == frag_val)) adjacent_positions.push_back(Position(x  ,y-1));

    return adjacent_positions;
  }


  // gets all adjacent pixels but NOT those diagonally connected
  vector<Position> get_all_adjacent_nodiag(size_t x,size_t y,uint32_t bg_val) {

    vector<Position> adjacent_positions;

    if(is_on_image(x+1,y  ) && (get(x+1,y  ) != bg_val)) adjacent_positions.push_back(Position(x+1,y  ));
    if(is_on_image(x-1,y  ) && (get(x-1,y  ) != bg_val)) adjacent_positions.push_back(Position(x-1,y  ));
    if(is_on_image(x  ,y+1) && (get(x  ,y+1) != bg_val)) adjacent_positions.push_back(Position(x  ,y+1));
    if(is_on_image(x  ,y-1) && (get(x  ,y-1) != bg_val)) adjacent_positions.push_back(Position(x  ,y-1));

    return adjacent_positions;
  }


  // gets all adjacent pixels but NOT those diagonally connected (BGVALS)
  vector<Position> get_all_adjacent_nodiag_bg(size_t x,size_t y,uint32_t bg_val) {

    vector<Position> adjacent_positions;

    if(is_on_image(x+1,y  ) && (get(x+1,y  ) == bg_val)) adjacent_positions.push_back(Position(x+1,y  ));
    if(is_on_image(x-1,y  ) && (get(x-1,y  ) == bg_val)) adjacent_positions.push_back(Position(x-1,y  ));
    if(is_on_image(x  ,y+1) && (get(x  ,y+1) == bg_val)) adjacent_positions.push_back(Position(x  ,y+1));
    if(is_on_image(x  ,y-1) && (get(x  ,y-1) == bg_val)) adjacent_positions.push_back(Position(x  ,y-1));

    return adjacent_positions;
  }


  void zero_image() {
    m_image_data = vector<uint32_t>(m_width*m_height,0);
  }

  Image get_region_image(size_t min_x,size_t max_x,size_t min_y,size_t max_y) {

    Image i;
    i.m_width  = max_x-min_x;
    i.m_height = max_y-min_y;

    i.zero_image();

    for(size_t x=min_x;x<max_x;x++) {
      for(size_t y=min_y;y<max_y;y++) {
        i.set(x-min_x,y-min_y,get(x,y));
      }
    }

    return i;
  }

  // Extracts the contiguous region representing this fragment.
  // This method operates in 2 parts:
  // 1. Find the region extreme values, cut this region out and creating a new image.
  // 2. Remove the region from the image
  Image grow_and_extract(size_t x,size_t y,int32_t bg_val) {
  

    //cout << "grow and extract: " << x << "," << y << endl;
    Image this_copy = *this;
    this_copy.zero_image();

    size_t min_x = 10000000000;
    size_t max_x = 0;
    size_t min_y = 10000000000;
    size_t max_y = 0;

    this_copy.set(x,y,get(x,y));
    set(x,y,bg_val);
    // Find region extremities
    vector<Position> positions = get_all_adjacent(x,y,bg_val);

    if(positions.size() == 0) {Image i; i.m_width=0; i.m_height=0; return i;} // It was a single pixel image, which breaks the code below, return a zero sized image...

    size_t i=0;
    for(;positions.size() != 0;) {
      vector<Position> o_pos;
      for(size_t n=0;n<positions.size();n++) {
	if(positions[n].m_x < min_x) min_x = positions[n].m_x;
	if(positions[n].m_x > max_x) max_x = positions[n].m_x;
	if(positions[n].m_y < min_y) min_y = positions[n].m_y;
	if(positions[n].m_y > max_y) max_y = positions[n].m_y;

        o_pos.push_back(Position(positions[n].m_x,positions[n].m_y));

        if(get(positions[n].m_x,positions[n].m_y) != bg_val) {
          this_copy.set(positions[n].m_x,positions[n].m_y,get(positions[n].m_x,positions[n].m_y));
        }
	set(positions[n].m_x,positions[n].m_y,bg_val);
      }

      positions.clear();

      for(size_t n=0;n<o_pos.size();n++) {
        vector<Position> p = get_all_adjacent(o_pos[n].m_x,o_pos[n].m_y,bg_val);
        
        // This ensures that the positions don't get used again, by setting them to bgval
        for(size_t i=0;i<p.size();i++) {
          if(get(p[i].m_x,p[i].m_y) != bg_val) {
            this_copy.set(p[i].m_x,p[i].m_y,get(p[i].m_x,p[i].m_y));
          }
          set(p[i].m_x,p[i].m_y,bg_val);
        }

        positions.insert(positions.begin(),p.begin(),p.end());
      }

      if((i%1000 == 0) && (i!=0)) cerr << "Growing region, current bounds: " << min_x << "," << max_x << " " << min_y << "," << max_y << endl;
      i++;
    }

    return this_copy.get_region_image(min_x,max_x,min_y,max_y);
  }

  Image trim(uint32_t bg_val) {
    
    size_t min_x = 10000000000;
    size_t max_x = 0;
    size_t min_y = 10000000000;
    size_t max_y = 0;

    bool nothing=true;
    for(size_t x=0;x<m_width;x++) {
      for(size_t y=0;y<m_height;y++) {
        if(get(x,y) != bg_val) {
          if(x<min_x) min_x = x;
          if(x>max_x) max_x = x;

          if(y<min_y) min_y = y;
          if(y>max_y) max_y = y;
          nothing = false;
        }
      }
    }

    if(nothing==true) {
      min_x = 0;
      max_x = 0;
      min_y = 0;
      max_y = 0;
    }

    return get_region_image(min_x,max_x,min_y,max_y);
  }

  /// Performs image segmentation, contiguous regions that are not "value" are placed in the same segment.
  vector<Image> segment_on_value(uint32_t value) {

    Image this_image = *this; // create a copy of ourselves, this segmentation will be destructive.

    vector<Image> all_fragments;

    for(size_t x=0;x<m_width;x++) {
      for(size_t y=0;y<m_height;y++) {

        if(((y*m_width)+x)%1000 == 0) cerr << "Current segmentation position: " << x << "," << y << endl;

        if(this_image.get(x,y) != value) {
          //cout << "DBG get value was non bg" << endl;
          Image result = this_image.grow_and_extract(x,y,value);
          all_fragments.push_back(result);
          //if(result.pixel_count() > 20) {
            cerr << "Identified fragment, size: " << result.pixel_count() << " total: " << all_fragments.size() << endl;
          //  result.save_tiff_rgb(string("lastfrag.tif"));
          //}
        }
      }
    }

    return all_fragments;
  }
  
  size_t pixel_count() {
    return m_image_data.size();
  }

  // Return true if the cluster is greater than size
  bool cluster_greater(size_t x,size_t y,uint32_t bg_val,int size) {

    Position p(x,y);

    vector<Position> pixels;
    pixels.push_back(p);

    //cout << "cluster greater: " << x << " " << y << endl;
    size_t old_pixels_size = 1000000000;
    for(;pixels.size() != old_pixels_size;) {

      vector<Position> all_res;
      for(size_t n=0;n<pixels.size();n++) {
        vector<Position> res = get_all_adjacent_nodiag(pixels[n].m_x,pixels[n].m_y,bg_val);
        all_res.insert(all_res.begin(),res.begin(),res.end());
      } 
      all_res.insert(all_res.begin(),pixels.begin(),pixels.end());
      //for(size_t n=0;n<all_res.size();n++) cout << "all_res: " << all_res[n].m_x << " " << all_res[n].m_y << endl;

      sort(all_res.begin(),all_res.end());
      //for(size_t n=0;n<all_res.size();n++) cout << "all_srt: " << all_res[n].m_x << " " << all_res[n].m_y << endl;
      vector<Position>::iterator iend = unique(all_res.begin(),all_res.end());

      vector<Position> res_unique;
      res_unique.insert(res_unique.begin(),all_res.begin(),iend);

      // for(vector<Position>::iterator i=all_res.begin();i<all_res.end();i++) res_unique.push_back(*i);
      old_pixels_size = pixels.size();
      pixels = res_unique;

      //for(size_t n=0;n<pixels.size();n++) cout << "pixel: " << pixels[n].m_x << " " << pixels[n].m_y << endl;
      //cout << "size: " << pixels.size() << endl;
      if(pixels.size() > size) return true;
    }

    return false;
  }


  // Return true if the cluster is greater than size
  bool frag_cluster_greater(size_t x,size_t y,uint32_t frag_val,uint32_t bg_val,int size) {

    Position p(x,y);

    vector<Position> pixels;
    pixels.push_back(p);

    //cout << "cluster greater: " << x << " " << y << endl;
    size_t old_pixels_size = 1000000000;
    for(;pixels.size() != old_pixels_size;) {

      vector<Position> all_res;
      for(size_t n=0;n<pixels.size();n++) {
        vector<Position> res = get_all_adjacent_nodiag_val(pixels[n].m_x,pixels[n].m_y,bg_val,frag_val);
        all_res.insert(all_res.begin(),res.begin(),res.end());
      }
      all_res.insert(all_res.begin(),pixels.begin(),pixels.end());
      //for(size_t n=0;n<all_res.size();n++) cout << "all_res: " << all_res[n].m_x << " " << all_res[n].m_y << endl;

      sort(all_res.begin(),all_res.end());
      //for(size_t n=0;n<all_res.size();n++) cout << "all_srt: " << all_res[n].m_x << " " << all_res[n].m_y << endl;
      vector<Position>::iterator iend = unique(all_res.begin(),all_res.end());

      vector<Position> res_unique;
      res_unique.insert(res_unique.begin(),all_res.begin(),iend);

      // for(vector<Position>::iterator i=all_res.begin();i<all_res.end();i++) res_unique.push_back(*i);
      old_pixels_size = pixels.size();
      pixels = res_unique;

      //for(size_t n=0;n<pixels.size();n++) cout << "pixel: " << pixels[n].m_x << " " << pixels[n].m_y << endl;
      //cout << "size: " << pixels.size() << endl;
      if(pixels.size() > size) return true;
    }

    return false;
  }


  void fill_hole(uint32_t bg_val) {
    for(int x=0;x<m_width;x++) {
      for(int y=0;y<m_height;y++) {

        if(get(x,y) == bg_val) {
          bool notisolated = frag_cluster_greater(x,y,bg_val,bg_val,20);

          if(notisolated == false) {
            //find a near pixel
            uint32_t adjcol=bg_val;
            int cx = x;
            int cy = y;
            for(int n=0;n<1000;n++) {
              cx += (rand()%3)-1;
              // cy += (rand()%3)-1; // this would be more accurate, but due to features of this image, only x may work better (test DARPA dataset)
              if(is_on_image(cx,cy)) adjcol = get(cx,cy);
              if(adjcol != bg_val) break;
            }
            set(x,y,adjcol);
          }
        }

      }
    }
  }

  void remove_fragments(uint32_t frag_col,uint32_t bg_val,int fragsize=30) {
    for(int x=0;x<m_width;x++) {
      for(int y=0;y<m_height;y++) {
        if(get(x,y) == frag_col) {
          vector<Position> all_p;

          bool notisolated = frag_cluster_greater(x,y,frag_col,bg_val,fragsize); // remove all pixel clusters < N pixels in size.

          if(notisolated == false) set(x,y,bg_val);
        }
      }
    }
  }

  void remove_isolated_pixels(uint32_t bg_val) {
    for(int x=0;x<m_width;x++) {
      for(int y=0;y<m_height;y++) {
        if(get(x,y) != bg_val) {
          vector<Position> all_p;

          bool notisolated = cluster_greater(x,y,bg_val,100); // remove all pixel clusters < 10 pixels in size.

          if(notisolated == false) set(x,y,bg_val);
        }
      }
    }
  }

  void remove_hangs(uint32_t bg_val) {
    for(int x=0;x<(((int)m_width)-1);x++) {
      for(int y=0;y<(((int)m_height)-1);y++) {
        uint32_t sq = get_square(x,y,bg_val);
        if((sq == 6) || (sq == 9)) {
          set(x  ,y,bg_val);
          set(x+1,y,bg_val);
          set(x  ,y+1,bg_val);
          set(x+1,y+1,bg_val);
        }
      }
    }
  }

  Image rotate(size_t o_x,size_t o_y,double radians,uint32_t bg_val) {
    Image out = rotate_source(o_x,o_y,radians,bg_val);
    out.fill_gaps(bg_val);
    return out;
  }

  Image rotate_source(size_t o_x,size_t o_y,double radians,uint32_t bg_val) {
 
    Image new_image; // image to copy data in to
 
    size_t pad;
    if(m_width > m_height) pad = m_width;
                      else pad = m_height;
    new_image.m_width  = pad*4;
    new_image.m_height = pad*4;
    new_image.zero_image();

    size_t x_min = 1000000000;
    size_t x_max = 0;
    size_t y_min = 1000000000;
    size_t y_max = 0;
 
    for(uint32_t cx=0;cx<m_width;cx++) {
      for(uint32_t cy=0;cy<m_height;cy++) {
        double r_x, r_y;

        if(get(cx,cy) != bg_val) {
          rotate_point(cx,cy,o_x,o_y,radians,r_x,r_y);

          size_t x_pos = r_x+(2*pad);
          size_t y_pos = r_y+(2*pad);
          new_image.set(x_pos,y_pos,get(cx,cy));

          if(x_pos > x_max) x_max = x_pos;
          if(x_pos < x_min) x_min = x_pos;
          if(y_pos > y_max) y_max = y_pos;
          if(y_pos < y_min) y_min = y_pos;
        }
      }
    }

    return new_image.get_region_image(x_min,x_max,y_min,y_max);
    //return new_image;
  }

  void fill_gaps(uint32_t bg_val) {
    for(size_t x=0;x<m_width;x++) {
      for(size_t y=0;y<m_height;y++) {
        if(get(x,y) == bg_val) {
         
          if(is_on_image(x-1,y-1) && is_on_image(x+1,y+1)) { 
            vector<uint32_t> adjs;
            adjs.push_back(get(x+1,y));
            adjs.push_back(get(x-1,y));
            adjs.push_back(get(x,y+1));
            adjs.push_back(get(x,y-1));

            int perform_set=0;
            for(size_t n=0;n<adjs.size();n++) if(adjs[n] == bg_val) perform_set++;

            if(perform_set < 2) {
              sort(adjs.begin(),adjs.end());
              set(x,y,adjs[2]);          
            }
          }
        }
      }
    }
  }

  Image rotate_dest(size_t o_x,size_t o_y,double radians,uint32_t bg_val) {

    Image new_image; // image to copy data in to

    int32_t pad;
    if(m_width > m_height) pad = m_width;
                      else pad = m_height;
    new_image.m_width  = pad*4;
    new_image.m_height = pad*4;
    new_image.zero_image();

    for(int32_t cx=0-pad;cx<pad;cx++) {
      for(int32_t cy=0-pad;cy<pad;cy++) {
        double r_x, r_y;
        rotate_point(cx,cy,o_x,o_y,0-radians,r_x,r_y);
        if(is_on_image(r_x,r_y)) new_image.set(cx+pad,cy+pad,get(r_x,r_y));
      }
    }

    return new_image;
  }

  bool is_line_concave(size_t y,uint32_t bg_val) {
    bool was_inside   = false;
    bool went_outside = false;
    bool line_concave = false;


    // Clear out any isolated pixels, they are probably artifacts TODO: THIS COULD BE A SEPERATE IMAGE CLEANING FUNCTION...

    for(size_t x=1;x<(m_width-1);x++) {
      int pix[3];
      pix[0] = get(x-1,y);
      pix[1] = get(x  ,y);
      pix[2] = get(x+1,y);
      for(size_t n=0;n<3;n++) if(pix[n] != bg_val) pix[n] = 1; else pix[n] = 0;

      if((pix[1] != pix[0]) && (pix[1] != pix[2])) { set(x,y,get(x-1,y));}
    }

    for(size_t x=0;x<m_width;x++) { if(get(x,y) == bg_val) cout << "0"; else cout << "1"; }
    cout << endl;

    for(size_t x=0;x<m_width;x++) {

      if(get(x,y) != bg_val) {
	if(was_inside && went_outside) {
	  line_concave = true;
	}
      }

      if(get(x,y) != bg_val) was_inside = true;
 
      if(get(x,y) == bg_val) {
        if(was_inside) { went_outside = true; }
      }
    }

    if(line_concave) cout << "concave" << endl; else cout << "not concave" << endl;
    return line_concave;
  }

  int top_is_concave(uint32_t bg_val) {

    cout << "CHECKING TOP" << endl;
    size_t concave_count = 0;

    for(size_t y=0;y<(m_height/8);y++) {
      bool line_concave = is_line_concave(y,bg_val);
      if(line_concave) { concave_count++; }
    }

    cout << "TOP concave_count: " << concave_count << endl;
    return concave_count;
  }

  int bottom_is_concave(uint32_t bg_val) {

    cout << "CHECKING BOTTOM" << endl;
    size_t concave_count = 0;

    for(size_t y=(m_height-1);y>=(m_height-(m_height/8));y--) {
      bool line_concave = is_line_concave(y,bg_val);
      if(line_concave) concave_count++;
    }

    cout << "BOTTOM concave_count: " << concave_count << endl;
    return concave_count;
  }

  Image orientate(Image::fragment_type f_type,uint32_t bg_val) {

    size_t min_width=100000000;
    double min_width_radians = 0;

    for(double r=0;r<(2*3.141);r+=0.01) {
      Image i2 = rotate(m_width/2,m_height/2,r,0);


      if(i2.m_width < min_width) {
	min_width = i2.m_width;
	min_width_radians = r;
      }
    }

    cout << "Orientated at radians: " << min_width_radians << endl;
    Image rotated = rotate(m_width/2,m_height/2,min_width_radians,0);

//    if(!(rotated.top_is_concave(bg_val) && rotated.bottom_is_convex(bg_val))) {

    int top_conv_count = rotated.top_is_concave   (bg_val);
    int bot_conv_count = rotated.bottom_is_concave(bg_val);

    if(bot_conv_count > top_conv_count) {
      cerr << "Fragment upside-down" << endl;
      rotated.rotate180();
    }

    return rotated;
  }

  void swap_pixel(size_t x1,size_t y1,
                  size_t x2,size_t y2) {

    uint32_t t = get(x1,y1);
    set(x1,y1,get(x2,y2));
    set(x2,y2,t);
  }

  // Inplace 180 degree rotation.
  void rotate180() {

    size_t h;
    if((m_height%2) == 1) {
      h = (m_height/2)-1;
      size_t ymid = ((m_height-1)/2);

      for(size_t x=0;x<(m_width/2);x++) {
        swap_pixel(x,ymid,(m_width-x-1),ymid);
      }

    } else h = (m_height/2)-1;
    for(size_t x=0;x<m_width;x++) {
      for(size_t y=0;y<=h;y++) {
        swap_pixel(x,y,m_width-x-1,m_height-y-1);
      }
    }
  }

  Image erode(uint32_t bg_val) {

    Image i = *this;

    for(size_t x=1;x<(m_width-1);x++) {
      for(size_t y=1;y<(m_height-1);y++) {
        if(get(x,y) != bg_val) {

          if(!((get(x+1,y  ) != bg_val) && (get(x-1,y  ) != bg_val) &&
               (get(x  ,y+1) != bg_val) && (get(x  ,y-1) != bg_val))) {
            i.set(x,y,bg_val);
          }
        }
      }
    }

    return i;
  }

  Image clean(uint32_t bg_val) {

    Image i = *this;

    // Remove all image edge pixels
    for(size_t x=0;x<m_width ;x++) { i.set(x,m_height-1,bg_val); }
    for(size_t x=0;x<m_width ;x++) { i.set(x,0         ,bg_val); }

    for(size_t y=0;y<m_height;y++) { i.set(m_width-1,y,bg_val); }
    for(size_t y=0;y<m_height;y++) { i.set(0        ,y,bg_val); }

    int erode_size = 10;
    if(m_width < m_height) erode_size = (m_width/2)-35;
    cout << "es: " << erode_size << endl;
    if(erode_size > 15 ) erode_size = 15;
    if(erode_size < -20) erode_size = 1;
    if(erode_size < -10) erode_size = 2;
    if(erode_size <   0) erode_size = 5;

    cout << "Using erode size: " << erode_size << endl;

    for(int n=0;n<erode_size;n++) i = i.erode(bg_val);

    // We assume all the pixels in the eroded image are foreground pixels, we'll keep these.

    map<uint32_t,int> pixels;
    size_t pixelsused=0;
    for(size_t n=0;n<m_image_data.size();n++) {
      if(i.m_image_data[n] != bg_val) {pixelsused++; pixels[m_image_data[n]]++;}
    }
    cerr << "Eroded image contained: " << pixelsused << endl;
    pixels = grow_pixel_set(pixels,5);

    Image i2 = *this;
    i2.not_set_pixels(pixels,bg_val,0);
    i2.remove_isolated_pixels(bg_val);
    i2.fill_gaps(bg_val);
    i2.fill_gaps(bg_val);
    i2.fill_gaps(bg_val);
    i2.fill_gaps(bg_val);
    i2.fill_gaps(bg_val);
    i2.fill_gaps(bg_val);
    i2.erode(bg_val);

    i2 = i2.trim(bg_val);
    return i2;
  }

  // 0.5 == 50%
  void threshold(double fraction,uint32_t bg_val,uint32_t low_val,uint32_t high_val) {

    uint32_t frac_val = (fraction*((double)(((double)256)*256*256*256)));
    cout << "fraction : " << fraction << endl;
    cout << "frac_val : " << frac_val << endl;
    cout << " low_val : " <<  low_val << endl;
    cout << "high_val : " << high_val << endl;
    for(size_t n=0;n<m_image_data.size();n++) {
      if(m_image_data[n] != bg_val) {
        if(m_image_data[n] < frac_val) { m_image_data[n] = low_val; }
                                  else { m_image_data[n] = high_val;}
      }
    }
  }

  void locate_top_left(size_t &x,size_t &y,uint32_t bg_val) {
  
    size_t bestscore = 100000000;
    size_t best_x = 0;
    size_t best_y = 0;

    for(size_t cx=0;cx<(m_width/2);cx++) {

      for(size_t cy=0;cy<(m_height/2);cy++) {

        if(get(cx,cy) != bg_val) {
          int score = cx + cy;
          if(score < bestscore) {
            bestscore = score;
            best_x = cx;
            best_y = cy;
          }
        }
      }
    }

    x = best_x;
    y = best_y;

    for(;get_square(x,y,bg_val) == 15;) {
      x -= 1;
      y -= 1;
    }
  }

  void locate_bottom_left(size_t &x,size_t &y,uint32_t bg_val,fragment_type frag_type=frag_type_arrow) {
  
    size_t bestscore = 100000000;
    size_t best_x = 0;
    size_t best_y = 0;

    for(size_t cx=0;cx<(m_width/2);cx++) {

      for(size_t cy=m_height-1;cy>(m_height/2);cy--) {

        if(get(cx,cy) != bg_val) {
          int score = cx + (m_height-cy);
          if(frag_type==frag_type_rectangle) score = cx      + (m_height-cy);
          if(frag_type==frag_type_arrow    ) score = (cx*10) + (m_height-cy);
          if(score < bestscore) {
            bestscore = score;
            best_x = cx;
            best_y = cy;
          }
        }
      }
    }

    x = best_x;
    y = best_y;

    for(;get_square(x,y,bg_val) == 15;) {
      x -= 1;
      y += 1;
    }
  }

  void locate_top_right(size_t &x,size_t &y,uint32_t bg_val) {
  
    size_t bestscore = 100000000;
    size_t best_x = 0;
    size_t best_y = 0;

    for(size_t cx=(m_width/2);cx<m_width;cx++) {

      for(size_t cy=0;cy<(m_height/2);cy++) {

        if(get(cx,cy) != bg_val) {
          int score = (m_width-cx) + cy;
          if(score < bestscore) {
            bestscore = score;
            best_x = cx;
            best_y = cy;
          }
        }
      }
    }

    x = best_x;
    y = best_y;

    for(;get_square(x,y,bg_val) == 15;) {
      x += 1;
      y -= 1;
    }
  }


  void locate_bottom_right(size_t &x,size_t &y,uint32_t bg_val,fragment_type frag_type=frag_type_arrow) {
  
    size_t bestscore = 100000000;
    size_t best_x = 0;
    size_t best_y = 0;

    for(size_t cx=(m_width/2);cx<m_width;cx++) {

      for(size_t cy=m_height-1;cy>(m_height/2);cy--) {

        if(get(cx,cy) != bg_val) {
          int score;
          if(frag_type==frag_type_rectangle) score = (m_width-cx)      + (m_height-cy);
          if(frag_type==frag_type_arrow    ) score = ((m_width-cx)*10) + (m_height-cy);
          if(score < bestscore) {
            bestscore = score;
            best_x = cx;
            best_y = cy;
          }
        }
      }
    }

    x = best_x;
    y = best_y;

    for(;get_square(x,y,bg_val) == 15;) {
      x += 1;
      y += 1;
    }

  }


  // if only one way in and out, we call this a "bridge", basically means it only has two neighbours
  bool single_pixel_bridge(size_t x,size_t y,uint32_t bg_val) {

    vector<Position> p = get_all_adjacent_nodiag(x,y,bg_val);
  
    if(p.size() <= 2) return true; else return false;
  }


  int c_mt_x[16];
  int c_mt_y[16];
  int a_mt_x[16];
  int a_mt_y[16];

  void init_mt() {
    c_mt_x[0]  =  0; c_mt_y[0]  =  0;   // 0000
    c_mt_x[1]  =  1; c_mt_y[1]  =  0;   // 0001
    c_mt_x[2]  =  0; c_mt_y[2]  =  1;   // 0010
    c_mt_x[3]  =  1; c_mt_y[3]  =  0;   // 0011
    c_mt_x[4]  =  0; c_mt_y[4]  = -1;   // 0100
    c_mt_x[5]  =  0; c_mt_y[5]  = -1;   // 0101
    c_mt_x[6]  =  0; c_mt_y[6]  = -1;   // 0110
    c_mt_x[7]  =  0; c_mt_y[7]  = -1;   // 0111
    c_mt_x[8]  = -1; c_mt_y[8]  =  0;   // 1000
    c_mt_x[9]  =  1; c_mt_y[9]  =  0;   // 1001
    c_mt_x[10] =  0; c_mt_y[10] =  1;   // 1010
    c_mt_x[11] =  1; c_mt_y[11] =  0;   // 1011
    c_mt_x[12] = -1; c_mt_y[12] =  0;   // 1100
    c_mt_x[13] = -1; c_mt_y[13] =  0;   // 1101
    c_mt_x[14] =  0; c_mt_y[14] =  1;   // 1110
    c_mt_x[15] =  1; c_mt_y[15] =  0;   // 1111

    a_mt_x[0]  =  0; a_mt_y[0]  =  0;   // 0000
    a_mt_x[1]  =  0; a_mt_y[1]  =  1;   // 0001
    a_mt_x[2]  = -1; a_mt_y[2]  =  0;   // 0010
    a_mt_x[3]  = -1; a_mt_y[3]  =  0;   // 0011
    a_mt_x[4]  =  1; a_mt_y[4]  =  0;   // 0100
    a_mt_x[5]  =  0; a_mt_y[5]  =  1;   // 0101
    a_mt_x[6]  =  1; a_mt_y[6]  =  0;   // 0110 ??????????????
    a_mt_x[7]  = -1; a_mt_y[7]  =  0;   // 0111
    a_mt_x[8]  =  0; a_mt_y[8]  = -1;   // 1000
    a_mt_x[9]  =  0; a_mt_y[9]  =  1;   // 1001
    a_mt_x[10] =  0; a_mt_y[10] = -1;   // 1010
    a_mt_x[11] =  0; a_mt_y[11] = -1;   // 1011
    a_mt_x[12] =  1; a_mt_y[12] =  0;   // 1100
    a_mt_x[13] =  0; a_mt_y[13] =  1;   // 1101
    a_mt_x[14] =  1; a_mt_y[14] =  0;   // 1110
    a_mt_x[15] = -1; a_mt_y[15] =  0;   // 1111
  }

  uint32_t get_square(size_t x,size_t y,uint32_t bg_val) {

    int bb0,bb1,bb2,bb3;

    if(is_on_image(x  ,y  ) && (get(x  ,y  ) != bg_val)) bb3 = 1; else bb3 = 0;
    if(is_on_image(x+1,y  ) && (get(x+1,y  ) != bg_val)) bb2 = 1; else bb2 = 0;
    if(is_on_image(x  ,y+1) && (get(x  ,y+1) != bg_val)) bb1 = 1; else bb1 = 0;
    if(is_on_image(x+1,y+1) && (get(x+1,y+1) != bg_val)) bb0 = 1; else bb0 = 0;

//    cout << "bb3: " << bb3 << endl;
//    cout << "bb2: " << bb2 << endl;
//    cout << "bb1: " << bb1 << endl;
//    cout << "bb0: " << bb0 << endl;
    uint32_t i = (bb3 << 3) | (bb2 << 2) | (bb1 << 1) | bb0;

    return i;
  }


  //TODO FIX!!!! NOT THE MAX AND THIS IS CORRECT CHANGE NAME
  uint32_t get_square_max(size_t x,size_t y,uint32_t fg_val) {
     vector<uint32_t> b;

    if(is_on_image(x  ,y  )) b.push_back(get(x  ,y  ));
    if(is_on_image(x+1,y  )) b.push_back(get(x+1,y  ));
    if(is_on_image(x  ,y+1)) b.push_back(get(x  ,y+1));
    if(is_on_image(x+1,y+1)) b.push_back(get(x+1,y+1));

    uint32_t max=0;
    for(size_t n=0;n<b.size();n++) {if(b[n] == fg_val) max = fg_val;}
    return max;
  }

  bool walk_pixels(size_t start_x,
                   size_t start_y,
                   size_t end_x,
                   size_t end_y,
                   uint32_t bg_val,
                   uint32_t feature_val,
                   bool clockwise,
                   vector<uint32_t> &walk_values,
                   vector<Position> &walk_pos
                   ) {

    cout << "end: " << end_x << " " << end_y << endl;

    init_mt(); // initialise matching squares table
    uint32_t last_square;
    uint32_t current_square = -1;

    size_t x = start_x;
    size_t y = start_y;

    for(int n=0;n<100000;n++) {
      //cout << "at " << x << " " << y << endl;

      last_square    = current_square;
      current_square = get_square(x,y,bg_val);
    //cout << "current_square: " << current_square << endl;

      if(!clockwise) {
        x += a_mt_x[current_square];
        y += a_mt_y[current_square];
      } else {
        x += c_mt_x[current_square];
        y += c_mt_y[current_square];
      }

      if((x == end_x) && (y == end_y)) return true;

      walk_values.push_back(get_square_max(x,y,feature_val));
      walk_pos   .push_back(Position(x,y));
    }

    return false;
    
  }

  void linear_repair(uint32_t bgval,uint32_t fgval,int dist=6) {
    linear_repair_fwd(bgval,fgval,dist);
    linear_repair_rev(bgval,fgval,dist);
  }

  void linear_repair_fwd(uint32_t bgval,uint32_t fgval,int dist) {
    for(size_t y=0;y<m_height;y++) {
      int fg_count=0;
      int last_fg=100000;
      for(size_t x=0;x<m_width;x++) {

        if(get(x,y) == fgval) {
          if(last_fg < dist) {
            for(size_t nx=x-last_fg;nx<x;nx++) set(nx,y,fgval);
            last_fg=0;
          }
        }

        if(get(x,y) == 0) {
          if(last_fg < dist) {
            for(size_t nx=x-last_fg;nx<x;nx++) if(get(nx,y) != 0) set(nx,y,fgval);
          }
        }

        if(get(x,y) == fgval) last_fg=0;

        last_fg++;
      }
    }
  }

  void linear_repair_rev(uint32_t bgval,uint32_t fgval,int dist) {
    for(size_t y=0;y<m_height;y++) {
      int fg_count=0;
      int last_fg=100000;
      for(int x=m_width-1;x>=0;x--) {

        if(get(x,y) == fgval) {
          if(last_fg < dist) {
            for(size_t nx=x+last_fg;nx>x;nx--) set(nx,y,fgval);
            last_fg=0;
          }
        }

        if(get(x,y) == 0) {
          if(last_fg < dist) {
            for(size_t nx=x+last_fg;nx>x;nx--) if(get(nx,y) != 0) set(nx,y,fgval);
          }
        }

        if(get(x,y) == fgval) last_fg=0;

        last_fg++;
      }
    }
  }

  int get_merge_score(Image &other,int ox,int oy,uint32_t bg_val,int block_size=5) {


    cout << "other.m_width: " << other.m_width << endl;
    cout << "      m_width: " <<       m_width << endl;

    cout << "other.m_height: " << other.m_height << endl;
    cout << "      m_height: " <<       m_height << endl;
    int overlaps=0;

    for(int cx=0;cx<other.m_width;cx+=block_size) {
      for(int cy=0;cy<other.m_height;cy+=block_size) {

        // Position on other image cx+ox, cy+oy
        if (is_on_image(cx+ox,cy+oy) &&
            (get(cx+ox,cy+oy) != bg_val) &&
            (other.get(cx,cy) != bg_val)) {
          if(other.get(cx,cy) == get(cx+ox,cy+oy)) overlaps+=2;
                                              else overlaps+=3;
        }
      }
    }

    int good_adj=0;
    for(int cx=0;cx<other.m_width;cx+=block_size) {
      for(int cy=0;cy<other.m_height;cy+=block_size) {
        if(other.get(cx,cy) != bg_val)
        if(get_all_adjacent_nodiag(cx+ox,cy+oy,bg_val).size() != 0) good_adj++;


        vector<Position> p = get_all_adjacent_nodiag(cx+ox,cy+oy,bg_val);
        vector<Position> p1 = other.get_all_adjacent_nodiag(cx,cy,bg_val);
        if((p.size() < 4)  || (p1.size() < 4))
        for(size_t n=0;n<p.size();n++) {
          if(get(p[n]) == other.get(cx,cy)) {
            good_adj+=2;
            if(get(p[n]) < 0x99999999) {cout << "THIS GOOD" << endl; good_adj+=30;}
          } else {
            if((get(p[n]) != bg_val) && (other.get(cx,cy) != bg_val)) {
              if(get(p[n]) < 0x99999999) {cout << "THIS BAD" << endl; good_adj-=10;}
              if(other.get(cx,cy) < 0x99999999) {cout << "THIS BAD" << endl; good_adj-=10;}
            }
          }
          
        }

      }
    }
    cout << "Overlaps: " << overlaps << endl;
    cout << "Good adj: " << good_adj << endl;

    return good_adj-overlaps;
  }

  void overlay(Image &other,int ox,int oy,uint32_t bg_val) {

    // paste other at this position
    for(int cx=0;cx<other.m_width;cx++) {
      for(int cy=0;cy<other.m_height;cy++) {
        if(is_on_image(cx+ox,cy+oy)) if(other.get(cx,cy) != bg_val) set(cx+ox,cy+oy,other.get(cx,cy));
      }
    }
  }

  // returns merge quality
  Image merge(Image &other,uint32_t bg_val,int &max_score,int &max_score_ox,int &max_score_oy,int block_size=5,int border_size=80) {

    // extend this image

    // apply X and Y offsets to "other".
    // search for X and Y offset that minimises:
    //  * adjacent pixels of non-similar colour
    //  * overlapping pixels

    max_score = 0;
    max_score_ox=0;
    max_score_oy=0;

    Image t = add_border(border_size,0);

    for(int ox=0;ox<160;ox++) {
      for(int oy=0;oy<200;oy++) {
        cout << "checking offset: " << ox << " " << oy << endl;
        int score = t.get_merge_score(other,ox,oy,bg_val,block_size);
        if(score > max_score) {
          max_score = score;
          max_score_ox = ox;
          max_score_oy = oy;
        }
      }
    }

    cout << "Merge score   : " << max_score << endl;
    cout << "Merge position: " << max_score_ox << " " << max_score_oy << endl;
    t.overlay(other,max_score_ox,max_score_oy,bg_val);

    return t;
  }


  // TODO: should use bg_val to set background value (don't zero image set to bg_val)
  Image add_border(size_t size,uint32_t bg_val) {
    Image nimage;
    nimage.m_width  = (*this).m_width  + size*2;
    nimage.m_height = (*this).m_height + size*2;
    cout << "NImage size: " << nimage.m_width << " " << nimage.m_height << endl;

    nimage.zero_image();

    for(size_t cx=0;cx<m_width;cx++) {
      for(size_t cy=0;cy<m_height;cy++) {
        nimage.set(cx+size,cy+size,get(cx,cy));
      }
    }

    return nimage;
  }

  void blockify() {

    int block_size = 5;

    for(size_t x=0;x<m_width;x+=block_size) {
      for(size_t y=0;y<m_height;y+=block_size) {

        // get all pixels in the block
        vector<uint32_t> px;
        for(size_t cx=0;cx<block_size;cx++) {
          for(size_t cy=0;cy<block_size;cy++) {
            if(is_on_image(x+cx,y+cy)) px.push_back(get(x+cx,y+cy));
          }
        }

        // find representative value
        sort(px.begin(),px.end());
        uint32_t rep=px[px.size()/2];

        // set all pixels in block to value
        for(size_t cx=0;cx<block_size;cx++) {
          for(size_t cy=0;cy<block_size;cy++) {
            if(is_on_image(x+cx,y+cy)) set(x+cx,y+cy,rep);
          }
        }
        

      }
    }

  }



  vector<uint32_t> m_image_data;
  size_t m_width;
  size_t m_height;
};



#endif
