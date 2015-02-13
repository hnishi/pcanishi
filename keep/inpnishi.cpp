#include"nlib.h"

/*class Inp_nishi{
private:
public:
	string filename;
	Inp_nishi( const char *inputname );
};*/

Inp_nishi::Inp_nishi(const char *inputname){
   filename = inputname;
   //cout<<"input file name: "<<filename<<endl;
}

string Inp_nishi::read( string param ){
   ifstream ifs( filename.c_str() ); // ifstream cannot handle type string (filename); it must be char*
   if(ifs.fail()){//error handling
      cerr<<"ERROR: cannot open file; "<<filename<<endl;
      cerr<<"program was ended by Inp_nishi::read( "<<param<<" ) in inpnishi.cpp \n";
      //return "error";
      exit(1);
   }
   string found_param, buf_string;
   while(!ifs.eof()){// to read binary
      ifs >> buf_string;
      if( param == buf_string ){
         //cout<<"!!! found !!!"<<endl;
         ifs >> found_param;
	 cout<<param<<" <- "<<found_param<<endl;
	 if( found_param == string("space") ) return " ";
	 return found_param;
      }
   }
   cout<<"cannot find the parameter ("<<param<<") in parameter.txt\n";
   return "nothing";
}

