const int cost_del = 1;
const int cost_ins = 1;
const int cost_sub = 1;

int hamming_distance(string &s1,string &s2) {

  int distance=0;

  for(int n=0;(n<s1.size()) && (n<s2.size());n++) {
    if(s1[n] != s2[n]) { distance++; }
                  else { if(s1[n] == '1') {distance-=200;}}
  }

  return distance;
}

string trim(string str) {

  string strout;

  bool saw_1 = false;
  for(size_t n=0;n<str.size();n++) {
    if(str[n] == '1') saw_1 = true;
    if(saw_1) strout += str[n];
  }

  saw_1 = false;
  string strout1;
  for(int n=strout.size()-1;n>=0;n--) {
    if(strout[n] == '1') saw_1 = true;
    if(saw_1) strout1 += strout[n];
  }

  string frw;
  for(int n=strout1.size()-1;n>=0;n--) {
    frw += strout1[n];
  }

  return frw;
}

int edit_distance_r(const string &s1,const string &s2) {
  int n1 = s1.size();
  int n2 = s2.size();

  int* p = new int[n2+1];
  int* q = new int[n2+1];
  int* r;

  p[0] = 0;
  for( int j = 1; j <= n2; ++j )
    p[j] = p[j-1] + cost_ins;

  for( int i = 1; i <= n1; ++i )
    {
      q[0] = p[0] + cost_del;
      for( unsigned int j = 1; j <= n2; ++j )
        {
          int d_del = p[j] + cost_del;
          int d_ins = q[j-1] + cost_ins;
          int d_sub = p[j-1];
//          d_sub += ( s1[i-1] == s2[j-1] ? 0 : cost_sub );
          if((s1[i-1] == s2[j-1]) && (s1[i-1] == '1')) d_sub -= 2;
                                                  else d_sub += 0;
          if((s1[i-1] != s2[j-1])) {
            d_sub += cost_sub;
            if(s1[i-1] == '1') d_sub += 2;
            if(s2[j-1] == '1') d_sub += 2;
          }

          q[j] = std::min( std::min( d_del, d_ins ), d_sub );
      }
      r = p;
      p = q;
      q = r;
    }

  int tmp = p[n2];
  delete[] p;
  delete[] q;

  return tmp;
}

int black_count(string &s) {

  int count=0;
  for(size_t n=0;n<s.size();n++) {
    if(s[n] == '1') count++;
  }

  return count;
}

int calc_distance(const string &s1_in,const string &s2_in,int offset,bool fast) {

  string s1;
  string s2;

  if(fast == true) {
    for(size_t n=0;n<s1_in.size();n++) { if((n%2)==0) s1 += s1_in[n]; }
    for(size_t n=0;n<s2_in.size();n++) { if((n%2)==0) s2 += s2_in[n]; }
    offset = offset/2;
  } else {
    s1 = s1_in;
    s2 = s2_in;
  }

  if(s1.size() < 2) return -1;
  if(s2.size() < 2) return -1;

  int s1_start;
  int s1_end;
  int s2_start;
  int s2_end;

  if(offset < 0) {
    s1_start = 0-offset;
    s1_end   = s1.size();

    s2_start = 0;
    s2_end   = s1.size()+offset;
  } else {
    s1_start = 0;
    s1_end   = s1.size() - offset;
  
    s2_start = offset;
    s2_end   = s2.size();
  }

  //cout << "s1 range: " << s1_start << " " << s1_end << endl;
  //cout << "s2 range: " << s2_start << " " << s2_end << endl;

  string s1_region = s1.substr(s1_start,s1_end-s1_start);
  string s2_region = s2.substr(s2_start,s2_end-s2_start);

  //int score = hamming_distance(s1_region,s2_region);
  int score = edit_distance_r(s1_region,s2_region);
  //cout << "offset: " << offset << endl;
  //cout << "score : " << score << endl;
  int overlap_size = 0;
  if(s1_region.size() < s2_region.size()) overlap_size = s1_region.size();
                                     else overlap_size = s2_region.size();

  //cout << "score: " << score << " overlap_size: " << overlap_size << endl;
  if(overlap_size < 10) return 0;

//  cout << "offset      : "  << offset << endl;  
//  cout << "score       : "  <<  score << endl;
//  cout << "score scaled: "  << (100-score);
//  cout << "over  scaled: "  <<  100*((double)overlap_size/(double)s1.size()) << endl;

  return (100-score) + 100*((double)overlap_size/(double)s1.size());// + black_count(s1_region) + black_count(s2_region);
}



// overlappy
int edit_distance(const string &s1,const string &s2) {

  //cout << "edit_distance here" << endl;
  int offset=(0-s1.size())+1;
  offset += 20; // at least 20px overlap

  int max_score = -1;
  int max_score_offset = offset;
  //cout << "offset: " << offset << endl;
  //cout << "s1size: " << s1.size() << endl;
  int s1size = s1.size();
  int s2size = s2.size();
  int ssize;
  if(s1size < s2size) ssize = s1size-1; else ssize = s2size-1;

  ssize -= 20; // at least 20px overlap

  for(;offset<ssize;offset+=10) {
    int score = calc_distance(s1,s2,offset,true);

    if(score > max_score) {
      max_score = score;
      max_score_offset = offset;
    }
  }


  int overlap_max = max_score_offset+20;
  if((overlap_max < (ssize+20)) && ((max_score_offset-20) > (offset=(0-s1.size())+1) ))
  for(offset=max_score_offset-20;offset<overlap_max;offset++) {
    int score = calc_distance(s1,s2,offset,false);

    if(score > max_score) {
      max_score = score;
      max_score_offset = offset;
    }
  }

  cout << "max_score_offset: " << max_score_offset << endl;
  cout << "max_score       : " << max_score << endl;
  return max_score;
}
