#ifndef bamCov_H_
#define bamCov_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <list>
#include <map>
#include <iomanip>
#include "./comm.hpp"
#include "./DataClass.hpp"
#include "./GeneCount.hpp"
#include <cstdlib>
#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/cram.h>
#include <htslib/kseq.h>
#include <stdio.h>
#include <unordered_map>
#include <thread>


using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)


std::string bam_hdr_get_sample_name(bam_hdr_t* hdr)
{
	//	if ( !hdr ) 		error("[E:%s:%d %s] [E:%s:%d %s] Failed to read the BAM header",__FILE__,__LINE__,__FUNCTION__,__FILE__, __LINE__, __FUNCTION__);

	char *ptext = strdup(hdr->text);  
	const char *p = ptext; 
	const char *q, *r;
	int32_t n = 0;
	std::string sm;
	while( ( q = strstr(p, "@RG" ) ) != 0 ) 
	{
		p = q + 3;
		r = q = 0;
		if ( ( q = strstr(p, "\tID:" ) ) != 0 ) q += 4;
		if ( ( r = strstr(p, "\tSM:" ) ) != 0 ) r += 4;
		if ( r && q ) 
		{
			char *u, *v;
			for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
			for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
			*u = *v = '\0';
			if ( sm.empty() )
				sm = r;
			else if ( sm.compare(r) != 0 ) 
			{
				cerr << "Multiple sample IDs are included in one BAM file "<<sm<<"\t"<<r<<endl;
				//error("[E:%s:%d %s] [E:%s:%d %s] Multiple sample IDs are included in one BAM file - %s, %s",__FILE__,__LINE__,__FUNCTION__, __FILE__, __LINE__, __FUNCTION__, sm.c_str(), r);
			}
		}
		else break;
		p = q > r ? q : r;
		++n;
	}


	if ( sm.empty() ) 
	{	
		cerr<<"Sample ID information cannot be found"<<endl;
		//		warning("[W:%s:%d %s] Sample ID information cannot be found",__FILE__,__LINE__,__FUNCTION__);
	}
	free(ptext);
	return sm;
}

void bamCov_help()
{
	cout<<""
		"\n"
		"\tUsage: PopGeneCNV  -i  InSortBam.list -g gene.gff -o  outPrefix\n"
		"\n"
		"\t\t-i    <str>     input SAM/BAM files List, must be sort bam\n"
		"\t\t-g    <str>     inPut gff/gtf file to stat only gene CDS Depth/Coverage\n"
		"\t\t-o    <str>     prefix of output file\n"
		"\n"
		"\t\t-r    <str>     inPut the ref.fa file for cram\n"
		"\t\t-b    <str>     list of the regions of which the coverage and mean of depth would be given\n"
		"\n"
		"\t\t-v    <float>   the min variance(S2) of CNV to filter cnv[0.1]\n"
		"\t\t-p    <float>   the min cnv polymorphic to filter cnv [0.1]\n"
		"\n"
		"\t\t-l    <int>     the region extension length,default [300]\n"
		"\t\t-t    <int>     the thread Num(CPU) for deal bam[3]\n"
		"\t\t-s              the bam with No bai method\n"
		"\n"
		"\t\t-h              show more details for help [hewm2008 v1.02]\n"
		"\n";
}



void More_bamCov_help()
{
	bamCov_help();
	cout<<""
		"\t\t-f    <str>     Identify the gff/gtf file feature for gene region [CDS]\n"
		"\t\t-q    <int>     the quality to filter reads, default [0]\n"
		"\t\t-d              No filter the PCR/optical duplicate read(DUP)\n"
		"\t\t-e              No filter the secondary alignment read(SECONDARY)\n"
		"\t\t-c              No filter the not passing quality controls read(QCFAIL)\n"
		"\n";
}




//#define DEFAULT_DEPTH 64000

/*
   static int read_bam(void *data, bam1_t *b)
   {
   aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
   int ret;
   while (1)
   {
   ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->header, b);
   if ( ret<0 ) break;
   if ( b->core.flag & aux->flags ) continue;
   if ( (int)b->core.qual < aux->min_mapQ ) continue;
   break;
   }
   return ret;
   }
   */

int bamCov_help01(int argc, char **argv , In3str1v * paraFA04   )
{
	if (argc <=1 ) {bamCov_help();return 0;}
	int file_count=0;
	bool AA=true;   bool BB=true;   bool CC=true;

	for(int i = 1; i < argc ; i++){
		if(argv[i][0] != '-')
		{
			cerr << "command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");


		if (flag  == "InFile"  ||  flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			paraFA04->InStr1=argv[i];
		}
		else if (flag  ==  "OutPut"|| flag  ==  "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "Ref"|| flag  ==  "r" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->cram_reference=argv[i];
		}

		else if (flag  ==  "Gff"  ||   flag  =="g"   || flag  =="Gtf")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
			igzstream INGFF (argv[i],ifstream::in);
			if (INGFF.fail())
			{
				cerr << "open File error: "<<argv[i]<<endl;
				return  0;
			}

			string tmp ;
			for (int A=1 ; A<168 && (!INGFF.eof())  ; A++ )
			{
				getline(INGFF,tmp);			
				if (tmp.length()<2)
				{
					continue ;
				}
				if (tmp[0] == '#')
				{
					continue ;
				}
				if (tmp.find("Parent")!=string::npos)
				{
					paraFA04->InInt2=1 ;
				}
				else if (tmp.find("transcript_id")!=string::npos)
				{
					paraFA04->InInt2=2 ;
				}
			}

			INGFF.close();

			if ((paraFA04->InInt2)==0)
			{
				cerr<<"Error:InPut gff/gtf file can not find [Parent] or [transcript_id] to judge GFF or GTF,please check the -g file :"<<argv[i]<<endl;
				return 0;
			}
		}
		else if (flag  ==  "Bed"  ||   flag  =="b" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count++;
			(paraFA04->List).push_back(A);
		}
		else if (flag  ==  "S2"  || flag  =="v")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->MinS2=atof(argv[i]);
		}
		else if (flag  ==  "LogPi"  || flag  =="p")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->LogPi=atof(argv[i]);
		}
		else if (flag  ==  "CPU"  || flag  =="t")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->CPU=atoi(argv[i]);
		}
		else if (flag  ==  "feature"|| flag  ==  "f" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->CDS=argv[i];
		}
		else if (flag  ==  "length"  || flag  =="l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->readoverLen=atoi(argv[i]);
		}
		else if (flag  ==  "MinQ"  || flag  =="q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag  ==  "NoBai"   ||  flag  =="s")
		{
			paraFA04->TF=false ;
		}
		else if (flag == "help"  || flag  =="h")
		{
			More_bamCov_help();return 0;
		}
		else if (flag == "d")
		{
			AA=false;
		}
		else if (flag == "e")
		{
			BB=false;
		}
		else if (flag == "c")
		{
			CC=false;
		}
		else
		{
			cerr << "UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if (  (paraFA04->InStr1).empty()  || (paraFA04->InStr3).empty() )
		//if (  (paraFA04->InStr1).empty()  || (paraFA04->InStr3).empty() || (paraFA04->cram_reference).empty() )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	if (  (paraFA04->InStr2).empty()  && file_count==0 )
	{
		cerr<< "lack argument for the must"<<endl ;
		return 0;
	}

	(paraFA04->InStr3)=add_Asuffix(paraFA04->InStr3);



	paraFA04->flags = BAM_FUNMAP ;
	if  (AA)
	{
		paraFA04->flags = (paraFA04->flags | BAM_FDUP);
	}
	if  (BB)
	{
		paraFA04->flags = (paraFA04->flags | BAM_FSECONDARY);
	}
	if  (CC)
	{
		paraFA04->flags = (paraFA04->flags | BAM_FQCFAIL);
	}


	return 1 ;
}





int CheckUpDate(string InBin )
{

	time_t t = time(0); 
	char tmp[32]={'\0'};
	strftime(tmp, sizeof(tmp), "%Y%m%d",localtime(&t)); 
	string today=tmp; string outDay="20230511";

	string InConfi;
	string  tmpStr;
	Rand6str(tmpStr);
	string OutPut;
	string group=getCmdResult("id | awk  '{print $2}' ");
	if (group.find("bc_",4)!=string::npos)
	{
		string user=getCmdResult("id | awk  '{print $1}' ");
		if (user.find("heweiming",4)==string::npos )
		{
			string outDay2="20230511";
			if ( (user.find("shaaiee",4)==string::npos)  ||  (today> outDay2 ) )
			{
				bamCov_help();
				cout <<"\t\t试用期结束 or 受权等其它限制,用户[ "<<user<<" ] 不可以用\n";
				cout<<"\t\tcontact  hewm2008 , $66.88 to get the latest version for 3 months\n\n";
				return 0;
			}
		}
	}
	else if (today > outDay )
	{
		if(access(InBin.c_str(), W_OK) == 0)
		{
			InConfi="wget https://raw.githubusercontent.com/hewm2008/PopGeneCNV/main/bin/PopGeneCNV  -O " +InBin ;
			cout <<"Download..."<<InConfi<<endl;
			tmpStr="/tmp/Rect.tmp"+tmpStr;
			InConfi="wget https://raw.githubusercontent.com/hewm2008/PopGeneCNV/main/bin/PopGeneCNV  -T 30 -t 1 -O  "+tmpStr ;
			int status = std::system(InConfi.c_str());
			//			cerr<<InConfi<<"\t"<<status<<endl;
			if (status < 0)
			{
				std::cout << "Error: " << strerror(status) << '\n';
				return 0;
			}
			else if (status ==0)
			{
				InConfi="mv  "+tmpStr+"  "+InBin+"; chmod  755 " + InBin ;
				std::system(InConfi.c_str());
			}
			else
			{
				InConfi="rm -rf "+tmpStr;
				std::system(InConfi.c_str());
				cout <<"\t\t下载超时,请手动到下面网址下载最新版本;\n";
				cout <<"\t\t由于程序强制定期升级 or 试用期结束 or 受权等其它限制，本程序已经作费, hewm2008\n";
				cout <<"\t\tDownload the latest version : https://github.com/hewm2008/PopGeneCNV or contact  hewm2008\n";
				return 0;
			}
		}
		else
		{
			cout <<"\t\t程序:\t"<<InBin<<"\t用户没有权限覆盖,请手动到下面网址下载最新版本;\n";
			cout <<"\t\t由于程序强制定期升级 or 试用期结束 or 受权等其它限制，本程序已经作费, hewm2008\n";
			cout <<"\t\tDownload the latest version : https://github.com/hewm2008/PopGeneCNV or contact  hewm2008\n";
			return 0;
		}
	}

	return 1;
}


void ProDealChrBambai (string  & BamPath , In3str1v * paraFA04   ,  map <int,map <int,int> >  & RegionMerger , unsigned short int **  depth, int  ChrNum  ,string &  SampleName)
{

	htsFile *fphts;
	sam_hdr_t *headerAA;
	int64_t rcnt;

	int n=1;
	int i=0;
	hts_idx_t *idx=NULL;
	int64_t *cnt ;
	int ret ;

	int min_mapQ= (paraFA04->InInt);
	fphts = hts_open(BamPath.c_str(), "r");

	if(fphts->format.format == htsExactFormat::cram)
	{
		if ((paraFA04->cram_reference).empty())
		{
			cerr<<"InPut cram file must add -r ref.fa as together\n";
			return ;
		}
		const char* ref_file = (paraFA04->cram_reference).c_str();
		hts_set_fai_filename(fphts, ref_file);
	}




	//cout <<"Begin Bam with bai :"<<BamPath<<"\t for the "<<ChrNum+1<<" Chr"<<endl;
	if (fphts)
	{
		//	idx=sam_index_load(fphts,BamPath.c_str());
		idx=sam_index_load3(fphts, BamPath.c_str(), NULL, HTS_IDX_SILENT_FAIL);
	}

	if (fphts == 0 || idx == 0 ) 
	{
		cerr<<"ERROR: failed to open bam or bam bai file "<<BamPath<<endl; 
		return ;
	}

	headerAA = sam_hdr_read(fphts);

	if (headerAA == NULL) { cerr<<"ERROR: failed to read header for bam file "<<BamPath<<endl; return  ;}

	SampleName=bam_hdr_get_sample_name(headerAA);
	//		cnt = calloc(n, sizeof(*cnt));
	//		n_plp = calloc(n, sizeof(int));
	//	plp = calloc(n, sizeof(bam_pileup1_t*));
	bam1_t *aln = bam_init1();
	map <int,map <int,int> > :: iterator  RegionIt =  RegionMerger.find(ChrNum) ;
	map <int,int> :: iterator MapSSEE =(RegionIt->second).begin() ;
	uint32_t flags = paraFA04->flags ;
	//(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
	uint32_t *cigar;
	hts_itr_t *iter;
	for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++ )
	{
		int tid=RegionIt->first;
		int64_t  beg=(MapSSEE->first);
		int64_t end=(MapSSEE->second);
		bam_mplp_t mplp;
		int	pos;
		iter = sam_itr_queryi(idx, tid, beg, end);
		rcnt = 0;

		while ((ret = sam_itr_next(fphts, iter, aln)) >= 0)
		{
			if ( aln->core.flag & flags )  { continue; }
			if ( (aln->core).qual < (paraFA04->InInt) )	{	continue ;	}

			cigar = bam_get_cigar(aln);
			int32_t StartRead=((aln->core).pos);
			for(int i=0; i < aln->core.n_cigar;++i)
			{
				int cig=bam_cigar_op(cigar[i]);
				int ncig = bam_cigar_oplen(cigar[i]);

				if  (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF )
				{
					int endTmp=StartRead+ncig;
					for (  ; StartRead<endTmp;StartRead++)
					{
						depth[(aln->core).tid][StartRead]++;
					}
					continue;
				}
				else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
				{
					StartRead=StartRead+ncig;
				}
			}
		}


	}

	//cout<<"Done the "<<ChrNum+1<<" Chr"<<endl;

	bam_destroy1(aln);
	if (idx) hts_idx_destroy(idx);
	if (iter) sam_itr_destroy(iter);
	//if (iter) hts_itr_destroy(iter);
	if (headerAA) sam_hdr_destroy(headerAA);


}


int SampleBamDeal (string  & BamPath , In3str1v * paraFA04,  map <int,map <int,int> >  & RegionMerger , unsigned short int ** depth, map <string,int> &  Chr2IntMap , bam_hdr_t *header , ogzstream &  OUTPP )
{
	//	bool QCFAIL=true;

	//#define BAM_CMATCH      0
	//#define BAM_CINS        1
	//#define BAM_CDEL        2
	//#define BAM_CREF_SKIP   3
	//#define BAM_CSOFT_CLIP  4
	//#define BAM_CHARD_CLIP  5
	//#define BAM_CPAD        6
	//#define BAM_CEQUAL      7
	//#define BAM_CDIFF       8

	uint32_t flags =  paraFA04->flags;
	//(BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);


	map <int,map <int,int> > :: iterator MergerIt;
	bool  *EndChr = new bool [(header->n_targets)];
	map <int,int> :: iterator  *ArryIt = new map <int,int> :: iterator [(header->n_targets)];
	map <int,int> :: iterator  *ArryItEnd = new map <int,int> :: iterator [(header->n_targets)];
	map <int,int> :: iterator MapSSEE ;
	map <int,map <int,int> > :: iterator  RegionIt ;

	for(int i = 0; i < (header->n_targets); i++)  
	{
		int CC=(header->target_len[i])+500;
		EndChr[i]=false ;
		MergerIt=RegionMerger.find(i);
		if (MergerIt==RegionMerger.end())
		{
			EndChr[i]=true ;
		}
		else
		{
			ArryIt[i]=(MergerIt->second).begin();
			ArryItEnd[i]=(MergerIt->second).end();
		}

		for (int32_t j =0 ; j< CC ; j++)
		{
			depth[i][j]=0;
		}
	}

	bam1_t *aln = bam_init1();
	uint32_t *cigar;




	string  SampleName;



	string bambai=BamPath+".bai";
	string crambai=BamPath+".crai";
	if ( ( ( access(bambai.c_str(), 0) == 0 )  ||  (access(crambai.c_str(), 0) == 0 )    && (paraFA04->TF ) )  )
		//	if ( ( access(bambai.c_str(), 0) == 0 )  && (paraFA04->TF ) )	
	{

		vector<thread> ThreadsVector ;
		RegionIt=RegionMerger.begin();
		while (RegionIt!= RegionMerger.end())
		{
			for (int i=0 ; i<(paraFA04->CPU) ;i++ )
			{
				if (RegionIt!= RegionMerger.end())
				{
					int  ChrNum=RegionIt->first;
					ThreadsVector.emplace_back( ProDealChrBambai , std::ref(BamPath), paraFA04, std::ref(RegionMerger) , depth , ChrNum , std::ref(SampleName));
					RegionIt++;
				}
			}
			int AA=ThreadsVector.size();
			for (int i=0 ; i<AA ; i++)
			{
				ThreadsVector[i].join();
			}
			ThreadsVector.clear();
		}


		for(  RegionIt=RegionMerger.begin() ; RegionIt!= RegionMerger.end(); RegionIt++ )
		{
			MapSSEE=(RegionIt->second).begin() ;

		}




	}


	else
	{




		bam_hdr_t *headerRR;
		samFile *BamInRR = hts_open(BamPath.c_str(), "r");
		if(BamInRR->format.format == htsExactFormat::cram)
		{
			if ((paraFA04->cram_reference).empty())
			{
				cerr<<"InPut cram file must add -r ref.fa as together\n";
				return 1;
			}
			const char* ref_file = (paraFA04->cram_reference).c_str();
			hts_set_fai_filename(BamInRR, ref_file);
		}


		headerRR = sam_hdr_read(BamInRR);

		SampleName=bam_hdr_get_sample_name(headerRR);

		cout <<"Begin file with No bai :"<<BamPath<<endl;
		bool ALLEnd=false;

		while (sam_read1(BamInRR, header, aln) >= 0)
		{
			if ( EndChr[(aln->core).tid])  {continue ;}	
			if ( (aln->core).qual < (paraFA04->InInt) )
			{
				continue ;
			}
			if ( aln->core.flag & flags ) continue;
			int32_t EndRead=bam_endpos(aln);
			if (EndRead <(ArryIt[(aln->core).tid]->first)) {continue;}

			int32_t StartRead=((aln->core).pos);
			if ( StartRead  > (ArryIt[(aln->core).tid]->second) )
			{

				for (ArryIt[(aln->core).tid]++  ; ArryIt[(aln->core).tid]!=ArryItEnd[(aln->core).tid] ;ArryIt[(aln->core).tid]++)
				{
					if ( StartRead <=  (ArryIt[(aln->core).tid]->second)  )
					{
						break;
					}
				}
				if  (ArryIt[(aln->core).tid]==ArryItEnd[(aln->core).tid])
				{
					EndChr[(aln->core).tid]=true;
					ALLEnd=true;
					for(int i = 0; i < (header->n_targets); i++)
					{
						if (!EndChr[i])
						{
							ALLEnd=false;
						}
					}
					if (ALLEnd)
					{
						break;
					}
				}
			}

			cigar = bam_get_cigar(aln);
			for(int i=0; i < aln->core.n_cigar;++i)
			{
				int cig=bam_cigar_op(cigar[i]);
				int ncig = bam_cigar_oplen(cigar[i]);
				if  (cig==BAM_CMATCH || cig==BAM_CEQUAL  || cig==BAM_CDIFF )
				{
					int endTmp=StartRead+ncig;
					for (  ; StartRead<endTmp;StartRead++)
					{
						depth[(aln->core).tid][StartRead]++;
					}
					continue;
				}
				else if  (cig==BAM_CDEL ||   cig==BAM_CREF_SKIP )
				{
					StartRead=StartRead+ncig;
				}
			}


		}

		sam_close(BamInRR);		
		bam_hdr_destroy(headerRR);

	}

	bam_destroy1(aln);

	cout <<"In.sort.Bam Read done : "<<BamPath<<endl;
	if (SampleName.empty())
	{
		vector<string> Tmp1;
		vector<string> Tmp2;
		split(BamPath,Tmp1,"\\");
		split(Tmp1[Tmp1.size()-1],Tmp2,".");
		SampleName=Tmp2[0];
	}

	//	llong Ascii[256] = {0};


	/*
	   map <int,string>  RefBase ;
	   gzFile fpRef;
	   kseq_t *seq;
	   int l;
	   fpRef = gzopen((paraFA04->cram_reference).c_str(), "r");
	   seq = kseq_init(fpRef);


	   while ((l = kseq_read(seq)) >= 0)
	   {
	//	stat_str_base(seq->seq.s , Ascii ,seq->seq.l );
	string chr=(seq->name.s);
	int ID=Chr2IntMap[chr];
	string seqBB=seq->seq.s;
	RefBase.insert( map <int,string>  :: value_type (ID,seqBB));
	}
	OUT.close();
	int GCGCArry[256] = {0};	GCGCArry['C']=1;  GCGCArry['c']=1;
	GCGCArry['G']=1;  GCGCArry['g']=1;
	*/
	//	llong GC_Aleng=Ascii['C']+Ascii['G'] +0 + Ascii['c']+Ascii['g'];
	//	llong eff_Acout=Ascii['A']+Ascii['a'] +0 + Ascii['T']+Ascii['g']+GC_Aleng;



	map <string,int> :: iterator  MapItChr2Int ;

	if  (!(paraFA04->List).empty())
	{
		string ListName=(paraFA04->List)[0];
		igzstream LIST (ListName.c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "open bed File error: "<<ListName<<endl;
			return  1;
		}

		//string StatBed=PrefixO+".bed.stat.gz";
		//ogzstream  OUTPP (StatBed.c_str());


		int Start ; int End ;  string ChrName ;
		ubit64_t   SS_Len =0;
		ubit64_t  SS_Cov=0;
		ubit64_t  SS_TotalD=0;

		//		OUTPP<<"#chr\tStart\tEnd\tCoverage%\tMeanDepth"<<endl;
		while(!LIST.eof())
		{
			string  line;
			getline(LIST,line);
			if (line.length()<=0)  {continue;}
			if (line[0] == '#')  { OUTPP<<line<<"\tSampleName\tCoverage%\tMeanDepth\n"; continue;}
			istringstream isone (line,istringstream::in);
			isone>>ChrName>>Start>>End;
			if  (Start > End )
			{
				cerr<<line<<"\tThis region may wrong\n"<<endl;
				continue;
			}
			Start--; End--;

			MapItChr2Int=Chr2IntMap.find(ChrName);

			if (MapItChr2Int==Chr2IntMap.end())
			{
				cerr<<line<<"\tThis region may wrong\n"<<endl;
			}
			else
			{
				int Length=End-Start+1;
				ubit64_t SumDepth=0;
				int NumCover=0;
				End++;
				for (int ii=Start; ii<End ; ii++)
				{
					if ( depth[MapItChr2Int->second][ii]>0 )
					{
						NumCover++;
						SumDepth+=(depth[MapItChr2Int->second][ii]);
					}
				}
				double Coverage=NumCover*100.0/Length ;
				double MeanDepth=SumDepth*1.0/Length ;
				OUTPP<<line<<"\t"<<SampleName<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<endl;
				SS_Cov+=NumCover;
				SS_Len+=Length;
				SS_TotalD+=SumDepth;
			}
		}
		LIST.close();
		double Coverage=SS_Cov*100.0/SS_Len ;
		double MeanDepth=SS_TotalD*1.0/SS_Len;
		OUTPP<<"##\t"<<SampleName<<"\tRegionLength:\t"<<SS_Len<<"\tCovered Site:\t"<<SS_Cov<<"\tCoverage%:\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth\t"<<MeanDepth<<endl;

		//		OUTPP.close();
	}



	if  ((paraFA04->InInt2)!=0)
	{

		//		string StatBed=PrefixO+".gene.stat.gz";
		//		ogzstream  OUTPP (StatBed.c_str());

		igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "open GFF/GTF File error: "<<(paraFA04->InStr1)<<endl;
			return  1;
		}

		int Start ; int End ;  string ChrName ;
		ubit64_t  SS_Len =0;
		ubit64_t  SS_Cov=0;
		ubit64_t  SS_TotalD=0;
		//		ubit64_t  SS_GCGC=0;
		unordered_map <string,int > GeneLength;
		unordered_map <string,int > GeneCover;
		unordered_map <string,int > GeneDepth;
		//		unordered_map <string,int > GeneGCGC;

		//		OUTPP<<"#chr\tStart\tEnd\tCoverage%\tMeanDepth"<<endl;

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
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
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


				Start--; End--;

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"\tThis region may wrong\n"<<endl;
				}
				else
				{
					int Length=End-Start+1;
					ubit64_t SumDepth=0;
					int NumCover=0;
					//					int GCGC=0;
					End++;
					for (int ii=Start; ii<End ; ii++)
					{
						if ( depth[MapItChr2Int->second][ii]>0 )
						{
							NumCover++;
							SumDepth+=(depth[MapItChr2Int->second][ii]);
						}
						//						GCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];
					}
					GeneLength[GeneID]+=Length;
					GeneCover[GeneID]+=NumCover;
					GeneDepth[GeneID]+=SumDepth;
					//					GeneGCGC[GeneID]+=GCGC;
					SS_Cov+=NumCover;
					SS_Len+=Length;
					SS_TotalD+=SumDepth;
					//					SS_GCGC+=GCGC;
				}
			}

		}
		else
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				line=replace_all(line,"\"","");
				line=replace_all(line,";","");
				istringstream isone (line,istringstream::in);
				string  flag , CDS , ZhengFu ,geneID ;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
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


				Start--; End--;

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"\tThis region may wrong\n"<<endl;
				}
				else
				{
					int Length=End-Start+1;
					ubit64_t SumDepth=0;
					int NumCover=0;
					//					int GCGC=0;
					End++;
					for (int ii=Start; ii<End ; ii++)
					{
						if ( depth[MapItChr2Int->second][ii]>0 )
						{
							NumCover++;
							SumDepth+=(depth[MapItChr2Int->second][ii]);
						}

						//						GCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];
					}
					GeneLength[GeneID]+=Length;
					GeneCover[GeneID]+=NumCover;
					GeneDepth[GeneID]+=SumDepth;
					//					GeneGCGC[GeneID]+=GCGC;
					SS_Cov+=NumCover;
					SS_Len+=Length;
					SS_TotalD+=SumDepth;
					//					SS_GCGC+=GCGC;
				}
			}
		}
		LIST.close();


		//double RefGCration=GC_Aleng*100.0/eff_Acout;

		//OUTPP<<"##RefGenome_EffLength:\t"<<eff_Acout<<"\tRefGenome_GCBase:\t"<<GC_Aleng<<setiosflags(ios::right)<<setprecision(2)<<"("<<RefGCration<<"%)"<<endl;


		/*
		   map <int ,double > GCDepth;
		   map <int, int > SumCount;
		   int GeneCount=0;
		   double GeneSumDepth=0;
		   for (auto iter = GeneLength.begin(); iter != GeneLength.end(); ++iter)
		   {
		   string GeneID=iter->first;
		   int Length=GeneLength[GeneID];
		   int SumDepth=	GeneDepth[GeneID];
		   int GCBase=	GeneGCGC[GeneID] ;
		   double MeanDepth=SumDepth*1.0/Length;
		   int  GeneGC=int(((GCBase*100)/Length));
		   SumCount[GeneGC]++;
		   GCDepth[GeneGC]+=MeanDepth;
		   GeneCount++;
		   GeneSumDepth+=MeanDepth;
		   }

		   double RDglobal=GeneSumDepth/GeneCount;
		   for (auto iter = GCDepth.begin(); iter !=GCDepth.end(); ++iter)
		   {
		   iter->second=iter->second/(SumCount[iter->first]);
		   }


		   OUTPP<<"#GeneID\tLength\tGCBase\tCoveredSite\tTotalDepth\tGC%\tCoverage%\tMeanDepth\tGeneCorRD\n";

*/

		OUTPP<<"#GeneID\tSampleName\tLength\tCoveredSite\tTotalDepth\tCoverage%\tMeanDepth\n";
		for (auto iter = GeneLength.begin(); iter != GeneLength.end(); ++iter)
		{
			string GeneID=iter->first;
			int Length=GeneLength[GeneID];
			int NumCover=GeneCover[GeneID];
			int SumDepth=	GeneDepth[GeneID];
			//	int GCBase=	GeneGCGC[GeneID] ;
			double Coverage=NumCover*100.0/Length ;
			double MeanDepth=SumDepth*1.0/Length;
			//	double GeneGC=GCBase*100.0/Length;
			//	int GCtmp=int(GeneGC);
			//	double RD_GC=GCDepth[GCtmp];
			//	double GeneCorDepth=RDglobal*MeanDepth/RD_GC;
			//	OUTPP<<GeneID<<"\t"<<Length<<"\t"<<GCBase<<"\t"<<NumCover<<"\t"<<SumDepth<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\t"<<GeneCorDepth<<endl;
			OUTPP<<GeneID<<"\t"<<SampleName<<"\t"<<Length<<"\t"<<NumCover<<"\t"<<SumDepth<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\t"<<Coverage<<"\t"<<MeanDepth<<endl;
		}


		double Coverage=SS_Cov*100.0/SS_Len ;
		double MeanDepth=SS_TotalD*1.0/SS_Len;
		//	double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

		OUTPP<<"##\t"<<SampleName<<"\tRegionLength:\t"<<SS_Len<<"\tCoveredSite:\t"<<SS_Cov<<"\tCoverage%:\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth\t"<<MeanDepth<<endl;
		//OUTPP<<"##RegionLength:\t"<<SS_Len<<"\tRegionGC\t"<<SS_GCGC<<"\tCoveredSite:\t"<<SS_Cov<<"\tCoverage%:\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth\t"<<MeanDepth<<"\tGC%\t"<<ALLGeneGCRation<<endl;
		//	OUTPP.close();
	}



	delete[] EndChr;
	delete[] ArryItEnd;




}



//int bamCov_main(int argc, char *argv[])
int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=-1;
	//paraFA04->TF=false ;
	if ((bamCov_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}
	string  InBin =argv[0];
	CheckUpDate(InBin);


	(paraFA04->InStr3)=(paraFA04->InStr3).substr(0,(paraFA04->InStr3).length()-3);
	string path=(paraFA04->InStr3);	
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	string PrefixO=path;

	if (ext == "stat" || ext == "bed")
	{
		PrefixO=path.substr(0,path.rfind('.'));
	}

	string  OUTStat=PrefixO+".gene.stat.gz";

	if  ((paraFA04->InInt2)==0)
	{
		OUTStat=PrefixO+".bed.stat.gz";
	}
	ogzstream  OUT (OUTStat.c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<OUTStat<<endl;
		delete  paraFA04 ; return  0;
	}
	vector <string>   BamFilevector ;
	ReadList ((paraFA04->InStr1)  ,BamFilevector  );

	string  BamPath=BamFilevector[0];
	if (BamPath.length()<=0)  
	{
		cerr<<"something wrong  at bam file"<<BamPath<<endl;
		return 1;
	}

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	if(BamIn->format.format == htsExactFormat::cram)
	{
		if ((paraFA04->cram_reference).empty())
		{
			cerr<<"InPut cram file must add -r ref.fa as together\n";
			return 1;
		}
		const char* ref_file = (paraFA04->cram_reference).c_str();
		hts_set_fai_filename(BamIn, ref_file);
	}

	header = sam_hdr_read(BamIn);

	map <string,int> Chr2IntMap; 
	for(int i = 0; i < (header->n_targets); i++) 
	{
		string ChrName=header->target_name[i];
		Chr2IntMap.insert( map <string,int>  :: value_type (ChrName,i));
	}

	sam_close(BamIn);		

	///////// swimming in the sky and flying in the sea ////////////



	map <int,map <int,int> > Region;
	map <int,map <int,int> > :: iterator  RegionIt ;
	map <int,int> :: iterator MapSSEE ;

	map <string,int> :: iterator  MapItChr2Int ;
	int Start ; int End ;  string ChrName ;
	int readoverLen=paraFA04->readoverLen;


	if  (!(paraFA04->List).empty())
	{
		string ListName=(paraFA04->List)[0];
		igzstream LIST (ListName.c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "open bed File error: "<<ListName<<endl;
			return  1;
		}


		while(!LIST.eof())
		{
			string  line;
			getline(LIST,line);
			if (line.length()<=0)  {continue;}
			if (line[0] == '#')  { continue;}
			istringstream isone (line,istringstream::in);
			isone>>ChrName>>Start>>End;
			if  (Start > End )
			{
				cerr<<line<<"\tThis region may wrong\n"<<endl;
				continue;
			}
			Start=Start-readoverLen;End=End+readoverLen;
			MapItChr2Int=Chr2IntMap.find(ChrName);

			if (MapItChr2Int==Chr2IntMap.end())
			{
				cerr<<line<<"\tThis region may wrong\n"<<endl;
			}
			else
			{
				RegionIt=Region.find(MapItChr2Int->second);
				if (RegionIt==Region.end())
				{
					map <int,int> Start2End ;
					Start2End[Start]=End;
					Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
				}
				else
				{
					MapSSEE=(RegionIt->second).find(Start);
					if (MapSSEE==(RegionIt->second).end())
					{
						(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
					}
					else
					{
						if  (End>(MapSSEE->second))
						{
							MapSSEE->second=End;
						}
					}
				}
			}
		}
		LIST.close();
	}









	if  ((paraFA04->InInt2)!=0)
	{

		string StatBed=PrefixO+".gene.stat.gz";
		ofstream  OUTPP (StatBed.c_str());

		igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "open GFF/GTF File error: "<<(paraFA04->InStr1)<<endl;
			return  1;
		}

		if ((paraFA04->InInt2)==1)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				istringstream isone (line,istringstream::in);
				string flag , CDS ;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End  ;
				Start=Start-readoverLen;End=End+readoverLen;
				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"\tThis region may wrong\n"<<endl;
				}
				else
				{

					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End ;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
						}
					}

				}
			}

		}
		else
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				istringstream isone (line,istringstream::in);
				string  flag , CDS ;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End  ;

				Start=Start-readoverLen;End=End+readoverLen;
				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"\tThis region may wrong\n"<<endl;
				}
				else
				{

					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
						}
					}

				}
			}
		}
		LIST.close();

	}



	map <int,map <int,int> > RegionMerger;

	map <int,map <int,int> > :: iterator MergerIt ;
	map <int,int> :: iterator  MergerMapSSEE  ;

	RegionIt=Region.begin() ;  	MapSSEE=(RegionIt->second).begin() ;



	for(  RegionIt=Region.begin() ; RegionIt!= Region.end(); RegionIt++ )
	{
		MapSSEE=(RegionIt->second).begin() ;
		int ChrInt = RegionIt->first ;
		Start=MapSSEE->first;
		End=MapSSEE->second;

		map <int,int> Start2End ;
		Start2End[Start]=End;
		RegionMerger.insert(map <int,map <int,int> > ::value_type(ChrInt,Start2End));
		MergerIt=RegionMerger.find(ChrInt);


		MapSSEE++;


		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++ )		
		{

			if (  (MapSSEE->first) >  End    )
			{
				Start=MapSSEE->first;
				End=MapSSEE->second;
				(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
			}
			else if  (  (MapSSEE->second) >  End    )
			{
				MergerMapSSEE=(MergerIt->second).end(); MergerMapSSEE--;
				MergerMapSSEE->second= End ;
				End=MapSSEE->second;
			}
		}
	}


	/*///


	  for(  RegionIt=RegionMerger.begin() ; RegionIt!= RegionMerger.end(); RegionIt++ )
	  {
	  MapSSEE=(RegionIt->second).begin() ;

	  for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++ )
	  {
	  cerr<<RegionIt->first<<"\t"<<(MapSSEE->first)<<"\t"<<(MapSSEE->second)<<"\n";
	  }
	  }



	  return 1;








	// *///





	///////// swimming in the sky and flying in the sea ////////////


	cout<<"begin new the memory ...\n";

	unsigned short int **depth = new unsigned short int*[(header->n_targets)];




	for(int i = 0; i < (header->n_targets); i++)  
	{
		int CC=(header->target_len[i])+500;
		depth[i] = new unsigned short int [CC];
	}



	cout<<"new the memory done"<<endl;
	int SampleNum=BamFilevector.size();
	for (int i=0 ;i < SampleNum ; i++)
	{
		BamPath=BamFilevector[i];
		SampleBamDeal(BamPath ,paraFA04, RegionMerger,depth, Chr2IntMap,header,OUT);
	}

	//�ͷſ��ٵ���Դ  
	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;

	bam_hdr_destroy(header);
	OUT.close();



	//	string OUTGeneHit=PrefixO+".GeneCNV.gz";

	GeneStat2Count(OUTStat ,PrefixO, paraFA04 ) ;



	////////////////   


	delete paraFA04 ;

	return 0;
}
#endif // bamCov_H_  //
///////// swimming in the sky and flying in the sea ////////////



