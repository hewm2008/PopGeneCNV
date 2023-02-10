#ifndef comm_H_
#define comm_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include  <zlib.h>
#include <htslib/kseq.h>
#include "gzstream/gzstream.c"
#include <sys/time.h>


/*///
#ifndef DBL_AMANT_ADIG
#define DBL_AMANT_ADIG 35
#endif

#define MAX_AFILE_ALIST_ALEN 512
#define MAX_ALIST_ANAME_ALEN 1024
#define MAX_ACHR_ANAME_ALEN 128
///*////
/*///
  const char decode[16] = {'N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N'};
  const int code[10] = {0,5,15,10,1,3,2,7,6,11};
  const int rev_Acode[16] = {0,4,6,5,4,1,8,7,6,8,3,9,5,7,9,2};
  const char abbv[17] = {'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};

//*/

using namespace std;

typedef long long llong ;
typedef unsigned long long ubit64_t;

KSEQ_INIT(gzFile, gzread)

	/*
	   template < class TT > 
	   string  Int2Str ( TT  A )
	   {
	   stringstream   sstrm ;
	   sstrm  <<  A ;
	   return  sstrm.str();
	   }
	   */

inline void  LogLackArg( string  flag )
{
	cerr << "\t\tLack Argument for [ -"<<flag<<" ]"<<endl;    
}

inline string add_Asuffix ( string path )
{
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "gz")
	{
		path=path+".gz" ; 
	}
	return path ;
}


string getCmdResult(const string &strCmd)
{
	char buf[10240] = {0};
	FILE *pf = NULL;

	if( (pf = popen(strCmd.c_str(), "r")) == NULL )
	{
		return "";
	}

	string strResult;
	while(fgets(buf, sizeof buf, pf))
	{
		strResult += buf;
	}

	pclose(pf);

	unsigned int iSize =  strResult.size();
	if(iSize > 0 && strResult[iSize - 1] == '\n')  // linux
	{
		strResult = strResult.substr(0, iSize - 1);
	}

	return strResult;
}





void Rand6str( string &  strResult )
{
	char chr[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
		'A', 'B', 'C', 'D', 'E', 'F', 'G',
		'H', 'I', 'J', 'K', 'L', 'M', 'N',
		'O', 'P', 'Q', 'R', 'S', 'T',
		'U', 'V', 'W', 'X', 'Y', 'Z',
		'a', 'b', 'c', 'd', 'e', 'f', 'g',
		'h', 'i', 'j', 'k', 'l', 'm', 'n',
		'o', 'p', 'q', 'r', 's', 't',
		'u', 'v', 'w', 'x', 'y', 'z'};

	struct timeval tv;
	gettimeofday(&tv, NULL);
	srand(tv.tv_sec + tv.tv_usec+((unsigned)time(NULL)));
	srand(time(NULL));
	for (int i=0; i<6; i++)
	{
		int idx = rand()%62;
		char buf=chr[idx];
		strResult=strResult+buf;
	}
}


string &  replace_all(string &  str,const  string &  old_Avalue,const string &  new_Avalue)
{
	while(true)   {
		string::size_type  pos(0);
		if(   (pos=str.find(old_Avalue))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}

string &   replace_all_distinct(string&   str,const   string &   old_Avalue,const   string &   new_Avalue)
{
	for(string::size_type   pos(0);   pos!=string::npos;   pos+=new_Avalue.length())   {
		if(   (pos=str.find(old_Avalue,pos))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}

inline string getID (string ID )
{
	string ext =ID.substr(0,ID.rfind('#')==string::npos ? ID.length() : ID.rfind('#')) ;
	return ext ;
}

void split(const string& str,vector<string>& tokens,  const string& delimiters = " ")
{
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}



int ReadList (string soaplist  ,   vector <string> & Soap_AStat  )
{
	igzstream LIST (soaplist.c_str(),ifstream::in); // igzstream
	int soapfilecout=0 ;
	if (!LIST.good())
	{
		cerr << "open List error: "<<soaplist<<endl;
		return  soapfilecout ;
	}
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0)  { continue  ; }
		Soap_AStat.push_back(line);
		soapfilecout++;
	}
	LIST.close();
	return  soapfilecout ;
}




string Int2Str (size_t A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

string Int2Str (int A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

string Int2Str (llong  A )
{
	stringstream   sstrm ;
	sstrm  <<  A ;
	return  sstrm.str();
}

inline void Swap ( int & x ,int & y)
{
	int tmp=y;
	y=x;
	x=tmp;
}


int stat_str_base(string  str , llong * Map , llong Leng )
{
	for(llong ix=0 ; ix<Leng ; ix++)
	{
		Map[str[ix]]++;
	}
	return 1 ;
}



#endif // comm_H_  ;

