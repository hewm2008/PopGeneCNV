
////////////////////////swimming in the sea & flying in the sky //////////////////


/*
 * DataClass.h
 *
 *  Created on: 2011-11-21
 *      Author: hewm@genomics.org.cn
 */

#ifndef DataClass_H_
#define DataClass_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "comm.hpp"

using namespace std;


///////////////// q_seq  for site ////////


class In3str1v {
	public:
		bool  TF ;
		int   InInt ;
		int   InInt2 ;
		uint32_t flags ;
		bool TF2 ;
		double InF ;
		int CPU ;
		int readoverLen ;
		string InStr1 ;
		string InStr2 ;
		string InStr3 ;
		string CDS ;
		string cram_reference ;
		vector <string> List ;
		double  LogPi ;
		double  MinS2 ;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			InStr3="";
			cram_reference="";
			TF=true ;
			TF2=true ;
			InInt=0 ;
			flags=0;
			InInt2=0;
			CPU=3;
			CDS="CDS";
			InF=0.0;
			readoverLen=300;
			LogPi=0.1;
			MinS2=0.1;
		}
};


////////swimming in the sky and flying in the sea *///////////





#endif /* DataClass_H_ */

//////////////// swimming in the sky and flying in the sea ////////////////



//////////////// swimming in the sky and flying in the sea ////////////////
