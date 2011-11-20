#include "tiffio.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "Image.h"
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <SDL/SDL.h>
#include "SDL_picofont.h"
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

  vector<string> get_top_left() {
    vector<string> all;
    sort(left_scores.begin(),left_scores.end());
    for(int n=left_scores.size()-1;n>(left_scores.size()-20);n--) all.push_back(left_scores[n].name);

    return all;
  }

  vector<string> get_top_right() {
    vector<string> all;
    sort(right_scores.begin(),right_scores.end());
    for(int n=right_scores.size()-1;n>(right_scores.size()-20);n--) all.push_back(right_scores[n].name);

    return all;
  }

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

  vector<Score> left_scores;
  vector<Score> right_scores;

  vector<Score> top_scores;
  vector<Score> bottom_scores;
};


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

///\brief reverses the endianness of a string
template <typename T> inline
void reverse_endian(T& t){
  unsigned char* res = reinterpret_cast<unsigned char*>(&t);
  unsigned char *temp = new unsigned char[sizeof(T)];
  for(int n=0;n<sizeof(T);n++) {temp[sizeof(T)-1-n] = res[n]; }
  for(int n=0;n<sizeof(T);n++) res[n] = temp[n];
  delete[] temp;
}

class QCView {
public:
  QCView(vector<string> filenames,string root,string adj_file_in) : m_filenames(filenames),m_adj_file(adj_file_in) {
    image_0_list = m_filenames;
    image_1_list = m_filenames;
    image_2_list = m_filenames;

    load_walks(root);

    m_adj_file_h = new ofstream(m_adj_file.c_str());
  }

  SDL_Surface *initialise_video() {
    if(SDL_Init(SDL_INIT_VIDEO)<0) {
      cerr << "Failed SDL_Init " << SDL_GetError() << endl;
      return NULL;
    }

    const SDL_VideoInfo *vidinfo = SDL_GetVideoInfo();

    int m_width  = vidinfo->current_w;
    int m_height = vidinfo->current_h-100;

    SDL_Surface *screen=SDL_SetVideoMode(m_width,m_height,32,SDL_ANYFORMAT);
    if(screen==NULL) {
      cerr << "Failed SDL_SetVideoMode: " << SDL_GetError() << endl;
      SDL_Quit();
      return NULL;
    }
    return screen;
  }

  void point(SDL_Surface *screen,int x,int y,int colour) {
    int bpp = screen->format->BytesPerPixel;
    Uint8 *p = (Uint8 *)screen->pixels + y * screen->pitch + x * bpp;
    if((x>=screen->w)||(y>=screen->h)||(x<0)||(y<0)) {}
    else {
      *(Uint32 *) p = colour;
    }
  }

  void render_image(size_t x,size_t y,Image i,uint32_t bg_val=0) {
    for(size_t cx=0;cx<i.m_width;cx++) {
      for(size_t cy=0;cy<i.m_height;cy++) {
        uint32_t p = i.get(cx,cy);

        //reverse_endian(p);  

        p = p & 0x00FFFFFF;

        //cout << p << endl;     

        if(p != bg_val) point(m_screen,cx+x,cy+y,p);
      }
    }
  }

  void render() {
    m_screen = initialise_video();

    current_image_1 = 0;
    current_image_2 = 0;

    i0_x=0;
    i0_y=0;

    i1_x=0;
    i1_y=0;
 
    i2_x=0;
    i2_y=0;

    i0_xd=0;
    i0_yd=0;

    i1_xd=0;
    i1_yd=0;

    i2_xd=0;
    i2_yd=0;

    i0.load_tiff(image_0_list[current_image_1]);
    i1.load_tiff(image_1_list[current_image_1]);
    i2.load_tiff(image_2_list[current_image_2]);

    for(size_t frame_number=0;;frame_number++) {

      SDL_Flip(m_screen);

      SDL_FillRect(m_screen,NULL,0);

      SDL_Event event;

      while(SDL_PollEvent(&event)) {
        process_key(event.key.keysym.sym,event.type);
      }

      render_image(i0_x,i0_y,i0);
      render_image(i1_x,i1_y,i1);
      render_image(i2_x,i2_y,i2);

      render_text(m_screen,image_0_list[current_image_0],0,0);
      render_text(m_screen,image_1_list[current_image_1],0,20);
      render_text(m_screen,image_2_list[current_image_2],0,40);

      i0_y += i0_yd;
      i0_x += i0_xd;
    
      i1_y += i1_yd;
      i1_x += i1_xd;

      i2_y += i2_yd;
      i2_x += i2_xd;
    }
    SDL_Quit();
  }

  void render_text(SDL_Surface *screen,
		   string text,
		   int x,
		   int y) {
    SDL_Color colour = {255, 255, 255};

    SDL_Rect destrec;

    // Render the text
    SDL_Surface *rendered_text = FNT_Render(text.c_str(), colour);

    destrec.x = x;
    destrec.y = y;
    destrec.w = rendered_text->w;
    destrec.h = rendered_text->h;

    if(x > screen->w) return;
    if(x < 0        ) return;
    if(y > screen->h) return;
    if(y < 0        ) return;

    if((x+destrec.w) > screen->w) return;
    if((x+destrec.w) < 0        ) return;
    if((y+destrec.h) > screen->h) return;
    if((y+destrec.h) < 0        ) return;


    SDL_BlitSurface(rendered_text, NULL, screen, &destrec);
  }

  void align_fragments() {
    i1.make_greyscale();
    i1.threshold(0.5,0,0x88888888,0xFFFFFFFF);
    i1.linear_repair(0,0x88888888);
    i1.remove_fragments(0x88888888,0,30);
    i1.blockify();
    i1.remove_fragments(0x88888888,0,110);

    i2.make_greyscale();
    i2.threshold(0.50,0,0x88888888,0xFFFFFFFF);
    i2.linear_repair(0,0x88888888);
    i2.remove_fragments(0x88888888,0,30);
    i2.blockify();
    i2.remove_fragments(0x88888888,0,110);


    int quality;
    int ox;
    int oy;
    Image nimage = i1.merge(i2,0,quality,ox,oy,5,80);
 
    i1_x=80;
    i1_y=80;
    i2_x=ox;
    i2_y=oy;
  }

 void find_best_neighbours() {
    // Find best scores for each fragment (currently left and right only)
    #pragma omp parallel for
    for(int j=0;j<fragments.size();j++) {
      int i = current_image_1;
      cout << "aligning: " << i << " " << j << endl;
      int lscore = edit_distance(fragments[i].leftside  ,fragments[j].rightside);
      int rscore = edit_distance(fragments[i].rightside ,fragments[j].leftside);
      int tscore = edit_distance(fragments[i].topside   ,fragments[j].bottomside);
      int bscore = edit_distance(fragments[i].bottomside,fragments[j].topside);


      #pragma omp critical
      {
        fragments[i].left_scores  .push_back(Score(lscore,j,fragments[j].filename));
        fragments[i].right_scores .push_back(Score(rscore,j,fragments[j].filename));
        fragments[i].top_scores   .push_back(Score(tscore,j,fragments[j].filename));
        fragments[i].bottom_scores.push_back(Score(bscore,j,fragments[j].filename));
      }
    }
 
    image_0_list = fragments[current_image_1].get_top_left();
    image_2_list = fragments[current_image_1].get_top_right();

    for(size_t n=0;n<image_0_list.size();n++) {find_and_replace(image_0_list[n],"edge_walks","orientated.tif");}
    for(size_t n=0;n<image_2_list.size();n++) {find_and_replace(image_2_list[n],"edge_walks","orientated.tif");}

    current_image_0 = 0;
    current_image_2 = 0;
    i0.load_tiff(image_0_list[current_image_0]);
    i2.load_tiff(image_2_list[current_image_2]);
  }

  void switch_to_clean() {

    string i0name = image_0_list[current_image_0];
    string i1name = image_1_list[current_image_1];
    string i2name = image_2_list[current_image_2];

    find_and_replace(i0name,"orientated","clean");
    find_and_replace(i1name,"orientated","clean");
    find_and_replace(i2name,"orientated","clean");

    i0.load_tiff(i0name);
    i1.load_tiff(i1name);
    i2.load_tiff(i2name);
    
  }

  void process_key(int key,int keystate) {
    int o_current_image_0 = current_image_0;  
    int o_current_image_1 = current_image_1;  
    int o_current_image_2 = current_image_2;  

    // key events

    if((key == SDLK_i) && (keystate == SDL_KEYDOWN)) { i0_y--; i0_yd=-1;}
    if((key == SDLK_k) && (keystate == SDL_KEYDOWN)) { i0_y++; i0_yd= 1;}
    if((key == SDLK_j) && (keystate == SDL_KEYDOWN)) { i0_x--; i0_xd=-1;}
    if((key == SDLK_l) && (keystate == SDL_KEYDOWN)) { i0_x++; i0_xd= 1;}


    if((key == SDLK_UP   ) && (keystate == SDL_KEYDOWN)) { i1_y--; i1_yd=-1;}
    if((key == SDLK_DOWN ) && (keystate == SDL_KEYDOWN)) { i1_y++; i1_yd= 1;}
    if((key == SDLK_LEFT ) && (keystate == SDL_KEYDOWN)) { i1_x--; i1_xd=-1;}
    if((key == SDLK_RIGHT) && (keystate == SDL_KEYDOWN)) { i1_x++; i1_xd= 1;}

    if((key == SDLK_w) && (keystate == SDL_KEYDOWN)) { i2_y--; i2_yd=-1;}
    if((key == SDLK_s) && (keystate == SDL_KEYDOWN)) { i2_y++; i2_yd= 1;}
    if((key == SDLK_a) && (keystate == SDL_KEYDOWN)) { i2_x--; i2_xd=-1;}
    if((key == SDLK_d) && (keystate == SDL_KEYDOWN)) { i2_x++; i2_xd= 1;}


    if((key == SDLK_r) && (keystate == SDL_KEYDOWN)) { align_fragments(); }
    if((key == SDLK_c) && (keystate == SDL_KEYDOWN)) { find_best_neighbours(); }
    if((key == SDLK_e) && (keystate == SDL_KEYDOWN)) { switch_to_clean(); }

    if((key == SDLK_n) && (keystate == SDL_KEYDOWN)) { (*m_adj_file_h) << image_0_list[current_image_0] << " " << image_1_list[current_image_1] << endl; }
    if((key == SDLK_m) && (keystate == SDL_KEYDOWN)) { (*m_adj_file_h) << image_1_list[current_image_1] << " " << image_2_list[current_image_2] << endl; }  


    if((key == SDLK_9)            && (keystate == SDL_KEYDOWN)) { current_image_0--; }
    if((key == SDLK_0)            && (keystate == SDL_KEYDOWN)) { current_image_0++; }
    if((key == SDLK_MINUS)        && (keystate == SDL_KEYDOWN)) { current_image_1--; }
    if((key == SDLK_EQUALS)       && (keystate == SDL_KEYDOWN)) { current_image_1++; }
    if((key == SDLK_LEFTBRACKET)  && (keystate == SDL_KEYDOWN)) { current_image_2--; }
    if((key == SDLK_RIGHTBRACKET) && (keystate == SDL_KEYDOWN)) { current_image_2++; }


    if(current_image_0 < 0) current_image_0 = 0;
    if(current_image_1 < 0) current_image_1 = 0;
    if(current_image_2 < 0) current_image_2 = 0;
    if(current_image_0 >= image_0_list.size()) current_image_0 = image_0_list.size()-1;
    if(current_image_1 >= image_1_list.size()) current_image_1 = image_1_list.size()-1;
    if(current_image_2 >= image_2_list.size()) current_image_2 = image_2_list.size()-1;
    if(o_current_image_0 != current_image_0) i0.load_tiff(image_0_list[current_image_0]);
    if(o_current_image_1 != current_image_1) i1.load_tiff(image_1_list[current_image_1]);
    if(o_current_image_2 != current_image_2) i2.load_tiff(image_2_list[current_image_2]);

    if(keystate == SDL_KEYUP) {cout <<"rset" << endl; i0_xd=0; i0_yd=0; i1_yd=0; i1_xd=0; i2_yd=0; i2_xd=0;}

    if((key == SDLK_q) && (keystate == SDL_KEYDOWN)) { SDL_Quit(); }

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

  void load_walks(string root) {
    // Load data
    vector<string> walk_filenames = get_all_edgewalks(root);

    for(int n=0;n<walk_filenames.size();n++) {
      cerr << "Loading: " << walk_filenames[n] << endl;

      string walk_filename = walk_filenames[n];

      fragments.push_back(FragmentInformation(walk_filename));
    }
  }

  SDL_Surface *m_screen;
  vector<string> image_0_list;
  vector<string> image_1_list;
  vector<string> image_2_list;


  vector<FragmentInformation> fragments;
  vector<string> m_filenames;
  ofstream *m_adj_file_h;
  string m_adj_file;

  int current_image_0; // left fragment
  int current_image_1; // centre fragment
  int current_image_2; // right fragment
  Image i0;
  Image i1;
  Image i2;

  int i0_x;
  int i0_y;

  int i1_x;
  int i1_y;

  int i2_x;
  int i2_y;

  int i0_xd;
  int i0_yd;

  int i1_xd;
  int i1_yd;

  int i2_xd;
  int i2_yd;
};

int main(int argc, char* argv[]) {

  cerr << "ds_qc <data folder> <adjacency file>" << endl;

  vector<string> image_filenames = get_all_orientated(string(argv[1]));

  QCView qcv(image_filenames,string(argv[1]),string(argv[2]));

  qcv.render();
}
