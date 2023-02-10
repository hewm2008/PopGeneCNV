#!/bin/sh
#$ -S /bin/sh
#Version1.0	hewm@genomics.org.cn	2023-01-16
echo Start Time : 
date
../bin/PopGeneCNV	-i	Bam.list	-g	/home/heweiming/new2022/KGBSA/01.Ref/MC.gff	-o	outPrefix -t 5	
echo End Time : 
date
