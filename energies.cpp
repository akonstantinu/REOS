//akonst02
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;


int main(){
    for(int i=0;i<1;i++){
      for(int j=0;j<20;j++){
	string en=std::to_string(double(i)+2.62958);
	string rot=std::to_string(0.999999-double(j)*0.02);
	system(("./nsss -f eosSA8 -e "+en+"e15 -r "+rot).c_str());
      }
    }
  return 0;
}
