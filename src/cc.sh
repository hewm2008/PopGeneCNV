#!/bin/sh
#$ -S /bin/sh
#Version1.0	 hewm@genomics.org.cn	 2022-12-14
echo Start Time : 
date
#/share/app/gcc-5.2.0/bin/g++	 --std=c++11	 -g	 	  PopGeneCNV.cpp	 -lcurl	 -lncurses	 -lhts	 -lm	 -lpthread	 -lboost_thread	 -lz	 -lrt	 	 -llzma	 -lbz2	 -I	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/include/	 -L	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/lib/	 -lcurl	 -lcrypto	 -ldl	 -L	 /home/heweiming/PS/dede/Lib/	 -o	 PopGeneCNV
/share/app/gcc-5.2.0/bin/g++	 --std=c++11	 -g	 -O3	  GeneCount.cpp	 -lcurl	 -lncurses	 -lhts	 -lm	 -lpthread	 -lboost_thread	 -lz	 -lrt	 -static	 -llzma	 -lbz2	 -I	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/include/	 -L	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/lib/	 -lcurl	 -lcrypto	 -ldl	 -L	 /home/heweiming/PS/dede/Lib/	 -o	 PopGeneCNV  -lpthread 
echo End Time : 
date
