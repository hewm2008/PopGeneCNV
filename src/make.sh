#!/bin/sh
#$ -S /bin/sh
#Version1.0	 hewm@genomics.org.cn	 2022-12-14
echo Start Time : 
date
echo /share/app/gcc-5.2.0/bin/g++	 --std=c++11	 -g	 -O3	  PopGeneCNV.cpp	 -lcurl	 -lncurses	 -lhts	 -lm	 -lpthread		 -lz	 	 -static	 -llzma	 -lbz2	 -I	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/include/	 -L	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/lib/	 -lcurl	 -lcrypto	 -ldl	 -L	 /home/heweiming/PS/dede/Lib/	 -o	 ../bin/PopGeneCNV  -lpthread 
/share/app/gcc-5.2.0/bin/g++	 --std=c++11	 -g	 -O3	  PopGeneCNV.cpp	 -lcurl	 -lncurses	 -lhts	 -lm	 -lpthread		 -lz	 	 -static	 -llzma	 -lbz2	 -I	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/include/	 -L	 /hwfssz4/BC_PUB/Software/08.Centos7/htslib-1.16/lib/	 -lcurl	 -lcrypto	 -ldl	 -L	 /home/heweiming/PS/dede/Lib/	 -o	 ../bin/PopGeneCNV  -lpthread 
echo End Time : 
date
