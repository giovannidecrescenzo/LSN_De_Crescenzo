#include <iostream>
#include <fstream>
#include <string>


using namespace std;

void read(string file){
   ifstream in;
   ofstream out;
   in.open("American_capitals.dat");
   out.open("American_coord.dat");
   while(!in.eof()){
      string city, state;
      double x, y;
      in >> city;
      in >> state;
      in >> x;
      in >> y;
      out << x << "\t" << y << endl;
   }
   in.close();
   out.close();
}

int main (int argc, char *argv[]){
   read("American_capitals.dat");
   return 0;
}