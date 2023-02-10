#ifndef GeneCount_H_
#define GeneCount_H_
#include "comm.hpp"
#include "DataClass.hpp"
#include <map>
#include <iomanip>
#include <unordered_map>



int CNV_Mat2VCF_Format(string  & INFile ,string & OutFile, In3str1v * paraFA04,string  & PASS)
{
	ogzstream OUT (OutFile.c_str());
	igzstream IN (INFile.c_str(),ifstream::in); // igzstream 
	if(!OUT.good())
	{
		cerr << "open OutFile error: "<<OutFile<<endl;
		exit(1);
	}
	if(!IN.good())
	{
		cerr << "open InputFile error: "<<INFile<<endl;
		exit(1);
	}


	unordered_map <string,int > GeneStart;
	unordered_map <string,int > GeneEnd;
	unordered_map <string,string > GeneChrID;
	unordered_map <string,int> :: iterator  UnMapItChr2Int ;

	igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in);
	string ChrName;

	if ((paraFA04->InInt2)==1)
	{
		while(!LIST.eof())
		{
			string  line;
			getline(LIST,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			istringstream isone (line,istringstream::in);
			string flag , CDS , ZhengFu ,geneID ;
			llong Start,End ;
			isone>>ChrName>>flag>>CDS ;
			if (CDS  != "CDS" )  { continue  ; }
			isone>>Start>>End>>flag>>ZhengFu>>flag>>geneID ;
			vector<string> inf;
			vector<string> Temp;
			split(geneID,inf,",;");
			split(inf[0],Temp,"=");
			string GeneID= Temp[Temp.size()-1] ;
			for ( int jk=1;  jk<inf.size() ; jk++)
			{
				vector<string> Temptmp2;
				split(inf[jk],Temptmp2,"=");
				if (Temptmp2[0] == "Parent")
				{
					GeneID= Temptmp2[Temptmp2.size()-1];
				}
			}
			UnMapItChr2Int=GeneStart.find(GeneID);
			if (UnMapItChr2Int==GeneStart.end())
			{
				GeneStart[GeneID]=Start;
				GeneChrID[GeneID]=ChrName;
			}
			else
			{
				if ( (UnMapItChr2Int->second) > Start )
				{
					(UnMapItChr2Int->second)= Start;
				}
			}

			UnMapItChr2Int=GeneEnd.find(GeneID);
			if (UnMapItChr2Int==GeneEnd.end())
			{
				GeneEnd[GeneID]=End;
			}
			else
			{
				if ( (UnMapItChr2Int->second) < End )
				{
					(UnMapItChr2Int->second)= End;
				}
			}

		}
	}
	else if ((paraFA04->InInt2)==2)
	{
		while(!LIST.eof())
		{
			string  line;
			getline(LIST,line);
			if (line.length()<=0 || line[0] == '#' )  { continue  ; }
			istringstream isone (line,istringstream::in);
			string flag , CDS , ZhengFu ,geneID ;
			llong Start,End ;
			isone>>ChrName>>flag>>CDS ;
			if (CDS  != "CDS" )  { continue  ; }

			isone>>Start>>End>>flag>>ZhengFu>>flag ;
			vector<string> inf;
			split(line,inf,"\t ");
			string GeneID=inf[9];
			for ( int jk=8;  jk<inf.size() ; jk+=2)
			{

				if (inf[0] == "transcript_id")
				{
					GeneID= inf[jk+1];
				}
			}

			UnMapItChr2Int=GeneStart.find(GeneID);
			if (UnMapItChr2Int==GeneStart.end())
			{
				GeneStart[GeneID]=Start;
				GeneChrID[GeneID]=ChrName;
			}
			else
			{
				if ( (UnMapItChr2Int->second) > Start )
				{
					(UnMapItChr2Int->second)= Start;
				}
			}

			UnMapItChr2Int=GeneEnd.find(GeneID);
			if (UnMapItChr2Int==GeneEnd.end())
			{
				GeneEnd[GeneID]=End;
			}
			else
			{
				if ( (UnMapItChr2Int->second) < End )
				{
					(UnMapItChr2Int->second)= End;
				}
			}

		}

	}

	LIST.close();

	
	
	OUT<<"##fileformat=VCF4.1\n"
		"##ALT=<ID=DEL,Description=\"Deletion\">\n"
		"##ALT=<ID=DUP,Description=\"Duplation\">\n"
		"##INFO=<ID=END,Number=1,Type=String,Description=\"End of the Gene Site\">\n"
		"##INFO=<ID=LogPi,Number=1,Type=String,Description=\"the LogPi of This Gene\">\n"
		"##INFO=<ID=S2,Number=1,Type=String,Description=\"variance of this CNV\">\n"
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

	string line ;
	getline(IN,line);
	vector<string> Temp;
	split(line,Temp," \t");
	for (int ii=3; ii<Temp.size(); ii++)
	{
		OUT<<"\t"<<Temp[ii];
	}
	OUT<<endl;



	
	
	// #GeneID	LogPi	S2	LB	LG	LH
	// MCh02g06373	2.38	0.13	0.00	0.76	0.00
	// MCh02g06524	0.17	0.11	0.37	1.02	0.28


//	string PASS=".";
	vector<string> inf;
	vector<double> infII;
	while(!IN.eof())
	{
		string line ;
		getline(IN,line);
		if (line.length()<=0) { continue ;}
		inf.clear();
		infII.clear();
		split(line,inf,"\t ");
		string ALT="<DEL>";
		int AAA=0; int BBB=0;
		infII.push_back(0);		infII.push_back(0);		infII.push_back(0);
		for (int ii=3; ii<Temp.size(); ii++)
		{
			double intCNV=atof(inf[ii].c_str());
			if  (intCNV<0.5) {AAA++;}
			else if  (intCNV>1.5) {BBB++;}
			infII.push_back(intCNV);
		}

		if (AAA<BBB)  {ALT="<DUP>";}

		OUT<<GeneChrID[inf[0]]<<"\t"<<GeneStart[inf[0]]<<"\t"<<inf[0]<<"\tN\t"<<ALT<<"\t.\t"<<PASS<<"\tEND="<<GeneEnd[inf[0]]<<";LogPi="<<inf[1]<<";S2="<<inf[2]<<"\tGT";

		for (int ii=3; ii<Temp.size(); ii++)
		{
			string GT="0/0";
			if  (infII[ii]>=0.5 && infII[ii]<=1.5 ) 
			{

			}
			else
			{GT="1/1";}
			OUT<<"\t"<<GT;
		}
		OUT<<"\n";
	}
// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Dad     S01     S02     S04     S05     S06     S07     S08     S12     S13     S14     S15
	OUT.close();
	IN.close();
}




int GeneStat2Count(string  & INFile ,string & OutFilePre, In3str1v * paraFA04)
{
	/*
	   int main(int argc,char *argv[])
	   {

	   if(argc!=3)
	   {
	   cout<<"\tVersion 1.0\thewm@genomics.org.cn\t2023-01-16\n";
	   cout<<argv[0]<<"\tInPut\tOutPut"<<endl;
	   return 1 ;
	   }
	   string INFile=argv[1];
	   string OutFilePre=argv[2];

	// */

	string OutFile=OutFilePre+".GeneCNV.Raw.mat.gz";
	string OutFile2=OutFilePre+".GeneDepth.mat.gz";
	string OutFile3=OutFilePre+".GeneCNV.Filter.mat.gz";


	ogzstream  OUT (OutFile.c_str());
	ogzstream  OUTDeth (OutFile2.c_str());
	ogzstream  OUT3 (OutFile3.c_str());
	igzstream IN (INFile.c_str(),ifstream::in); // igzstream 
	if(!OUT.good())
	{
		cerr << "open OutFile error: "<<OutFile<<endl;
		exit(1);
	}
	if(!IN.good())
	{
		cerr << "open InputFile error: "<<INFile<<endl;
		exit(1);
	}


	string GeneID,SampleName;
	int  Length,CoveredSite,TotalDepth;
	double	Coverage,MeanDepth;
	map <string,map <string,double> > :: iterator  ItMap  ;
	map <string,map <string,double> > GeneMap  ;
	map <string,bool> GeneIDMap;
	map <string,double> SNMap;
	map <string,double> :: iterator  FFF;
	map <string,bool> :: iterator EEE;

	////////////////////////swimming in the sea & flying in the sky //////////////////   

	while(!IN.eof())
	{
		string line ;
		getline(IN,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>GeneID>>SampleName>>Length>>CoveredSite>>TotalDepth>>Coverage>>MeanDepth;
		if ( GeneID[0] == '#' )
		{
			if( GeneID[1] == '#' )
			{
				vector<string> Temp;
				split(line,Temp," \t");
				MeanDepth=atof((Temp[9]).c_str());
				SNMap.insert(map <string,double> ::value_type(SampleName,MeanDepth));
			}
			continue;
		}

		EEE=GeneIDMap.find(GeneID);
		if (EEE==GeneIDMap.end())
		{
			GeneIDMap.insert(map <string,bool> ::value_type(GeneID,true));
		}

		ItMap=GeneMap.find(GeneID);
		if (ItMap==GeneMap.end())
		{
			map <string,double> Tmp;
			Tmp[SampleName]=MeanDepth;
			GeneMap.insert(map <string,map <string,double> > ::value_type(GeneID,Tmp));
		}
		else
		{
			(ItMap->second).insert(map <string,double> ::value_type(SampleName,MeanDepth));
		}

	}
	IN.close();
	//return 1;
	OUT<<"#GeneID\tLogPi\tS2";
	OUT3<<"#GeneID\tLogPi\tS2";
	OUTDeth<<"#GeneID";
	int SampleNumber=0;
	for (FFF = SNMap.begin(); FFF != SNMap.end(); ++FFF)
	{
		OUT<<"\t"<<FFF->first;
		OUT3<<"\t"<<FFF->first;
		OUTDeth<<"\t"<<FFF->first;
		SampleNumber++;
	}
	OUT<<endl;
	OUT3<<endl;
	OUTDeth<<endl;


	double *HitMean= new double [SampleNumber];
	double *PigLog= new double [SampleNumber];

	map <string,double> :: iterator  ItCCC  ;

	for ( EEE=GeneIDMap.begin();EEE!=GeneIDMap.end(); ++EEE)
	{
		OUT<<EEE->first<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2);
		OUTDeth<<EEE->first<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2);
		ItMap=GeneMap.find(EEE->first);
		int CountSN=0;
		double SumHit=0;
		for (FFF = SNMap.begin(); FFF != SNMap.end(); ++FFF)
		{
			ItCCC=(ItMap->second).find(FFF->first);
			double MeanDepth=ItCCC->second;
			double  MeanHit=MeanDepth/(FFF->second);
			//	OUT<<"\t"<<MeanHit;
			OUTDeth<<"\t"<<MeanDepth;
			HitMean[CountSN]=MeanHit;
			PigLog[CountSN]=log10(MeanHit*100+1);
			SumHit+=MeanHit;
			CountSN++;
		}
		SumHit=SumHit/SampleNumber;

		double S2=0;
		double S2Log=0;

		for (int ii=0 ;ii<SampleNumber ;  ii++)
		{
			S2+=(HitMean[ii]-SumHit)*(HitMean[ii]-SumHit);
			for (int jj=ii+1 ;jj<SampleNumber ;  jj++)
			{
				S2Log+=(PigLog[ii]-PigLog[jj])*(PigLog[ii]-PigLog[jj]);
			}
		}

		S2=S2/SampleNumber;
		S2Log=S2Log*2/(SampleNumber*(SampleNumber-1));
		OUT<<"\t"<<S2Log<<"\t"<<S2;

		for (int ii=0 ;ii<SampleNumber ; ii++)
		{
			OUT<<"\t"<<HitMean[ii];	
		}
		OUTDeth<<endl;
		OUT<<endl;


		//		if (S2Log>0.1  && S2>0.1)
		if (S2Log>(paraFA04->LogPi)  && S2>(paraFA04->MinS2))
		{			
			OUT3<<EEE->first<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<S2Log<<"\t"<<S2;
			for (int ii=0 ;ii<SampleNumber ; ii++)
			{
				OUT3<<"\t"<<HitMean[ii];	
			}
			OUT3<<endl;
		}
	}


	OUTDeth.close();	
	OUT.close();
	OUT3.close();	
	delete [] HitMean ;
	delete [] PigLog;



	string OutFile4=OutFilePre+".GeneCNV.Filter.vcf.gz";
	string PASS=".";
	CNV_Mat2VCF_Format( OutFile3 ,  OutFile4, paraFA04,PASS);
	OutFile4=OutFilePre+".GeneCNV.Raw.vcf.gz";
	PASS="PASS";
	CNV_Mat2VCF_Format( OutFile ,  OutFile4, paraFA04,PASS);



	return 0 ;
}

/*
int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=-1;
	string OutFile4=argv[2];
	string OutFile3=argv[1];
	string PASS=".";
	CNV_Mat2VCF_Format( OutFile3 ,  OutFile4, paraFA04,PASS);

}
*/

#endif // GeneCount_H_  ;
////////////////////////swimming in the sea & flying in the sky //////////////////


