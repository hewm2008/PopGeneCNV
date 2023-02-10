
# PopGeneCNV
A new simple and efficient software to Rapid detection of population gene CNV based on bam


###  1 Introduction
<b>PopGeneCNV</b>  ����bam (sam/cram)�ļ������ټ��Ⱥ�壨һ����Ʒ������Ʒ�����ڶ�̬��CNV�� ��ѡ�����������������Ʒ��bam�ļ���������Ϣ�ļ�(GFF/GTF���������Եõ� 1  Ⱥ��Ļ��������  2 Ⱥ��Ļ���CNV����  3 ��ѡ��Ⱥ�����CNV��̬�ԵĻ���
</br>
</br> keyword  : PopCNV ;&nbsp;&nbsp;&nbsp;&nbsp;    CNV  ;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; GFF/GTF ;&nbsp;&nbsp;&nbsp;&nbsp ;&nbsp;&nbsp;&nbsp;&nbsp;GeneCNV

</br>��Ҫ���ܣ�
</br>1  �Զ�����Ƿ��������bam.bai�ļ������ݻ������궨λbam�ļ�λ�ã����Լ�����߶�bam���ٶȣ���Ϊ���������ǻ������1%��ֻ��ȡ��Щ�����bam�Ἣ�����ٶȡ�   ͬʱ  ���� �Զ�ʶ��bam header��Ʒ��
</br>2  ���̣߳��ڴ���bam�����ļ�ʱ��Ĭ��3�����̲������ٶ�bam
</br>3  �Զ�ʶ��GFF����GTF,  ���ǵ���Щ���ֻ����intron���󣬳�������������޳�intron���򡣼���Ҫʶ������CDS���� 
</br>4  ������CNV��̬��ָ�꣬���ٹ��˳�Ⱥ���д���CNV����ĺ�ѡgene������
</br>5  �÷��򵥣����bam�ļ��б��ͻ�����Ϣ�ļ�(GFF/GTF) ,������һ����λ������� 
</br>6  ����и���Ʒ�Ļ�����Ⱥ͸��Ƕ��ļ����ͻ���CNV�����ļ����ʹ���Ⱥ���̬��CNV��ѡgene�ļ� 
</br>7  ��̬���룬linuxֱ�ӽ�ѹ����
</br>	


</br>�����Ǹ�һЩ�л��������������õģ�����С�׿����������ˡ�
</br>
</br><b>PopGeneCNV</b>  is based on the bam (sam/cram) file, the rapid detection population (one sample or multiple samples) exists in the waiting selection basis of CNV of various properties. That is, the bam file and gene file of the input sample and the FF/G file can be obtained (TF, G) to obtain the basic performance matrix of the 1 population, the basic CNV matrix of the 2 population, and the basis of 3 candidate populations existing in the CNV polymorphism



###  2 Download and Install
------------
The <b>new version</b> will be updated and maintained in <b>[hewm2008/PopGeneCNV](https://github.com/hewm2008/PopGeneCNV)</b>, please click below website to download the latest version
</br><p align="center"><b>[hewm2008/PopGeneCNV](https://github.com/hewm2008/PopGeneCNV)</b></p>

<b> 2.1. linux/MaxOS&nbsp;&nbsp;&nbsp;   [Download](https://github.com/hewm2008/PopGeneCNV/archive/v1.02beta.tar.gz)</b>
  
  </br> <b>2.2 Pre-install</b>
  </br> PopGeneCNV is for Linux/Unix/macOS only.
  </br>Before installing,please make sure the following pre-requirements are ready to use.
  </br> 1) g++   : g++ with [--std=c++11](https://gcc.gnu.org/) > 4.8+ is recommended
  </br> 2) htslib  : [htslib](https://github.com/samtools/htslib) > 1.15 is recommended
  </br> 3) lpthread    : [pthreads](https://man7.org/linux/man-pages/man7/pthreads.7.html) -lpthread is recommended

</br> <b>2.3 Install</b>
</br> Users can install it with the following options:
</br> Option 1: 
<pre>
        git clone https://github.com/hewm2008/PopGeneCNV.git
        cd PopGeneCNV;	chmod 755 -R bin/*
        ./bin/PopGeneCNV  -h 
</pre>


###  3 Parameter description
------------
</br><b>3.1 PopGeneCNV</b>
</br><b>3.1.1 Main parameter</b>

```php

./bin/PopGeneCNV

	Usage: PopGeneCNV  -i  InSortBam.list -g gene.gff -o  outPrefix

		-i    <str>     input SAM/BAM files List, must be sort sort.bam
		-g    <str>     inPut gff/gtf file to stat only gene CDS Depth/Coverage
		-o    <str>     prefix of output file
		
		-r    <str>     inPut the ref.fa file for cram
		-b    <str>     list of the regions of which the coverage and mean of depth would be given

		-v    <float>   the min variance(S2) of CNV to filter cnv[0.1]
		-p    <float>   the min cnv polymorphic to filter cnv [0.1]

		-q    <int>     the quality to filter reads, default [0]
		-l    <int>     the region extension length,default [300]
		-t    <int>     the thread Num(CPU) for deal bam[3]
		-s              the bam with No bai method

		-h              show more details for help [hewm2008 v1.02]


```
</br> brief description for function:
<pre>
	   #   �÷�һ���������������Ϊ ���������������� 

           ../bin/PopGeneCNV	-i	Bam.list	-g	~/01.Ref/MC.gff	-ooutPrefix -t 5


</pre>




</br><b>3.3 Output files</b>


|Module |    outFlie    |       Description                                                |
|:-----:|:-------------------------|:------------------------------------------------------------|
| List  |                          |                                                             |
|       |out.gene.stat.gz          |������Ʒ�����������Ⱥ͸��Ƕȵ�ͳ����Ϣ                             |
|       |out.GeneDepth.mat.gz      |������Ʒ�������ȣ������������(readcount��                                        |
|       |Out.GeneCNV.Raw.mat.gz    |��׼����Ⱥ�����CNV�ľ���                                      |
|       |Out.GeneCNV.Raw.mat.gz    |��׼����Ⱥ�����CNV��VCF��ʽ                                     |
|       |Out.GeneCNV.filter.mat.gz |��׼���󣬴���Ⱥ���̬CNV�ĺ�ѡ����                            |
|       |Out.GeneCNV.filter.vcf.gz |��׼���󣬴���Ⱥ���̬CNV�ĺ�ѡ����VCF��ʽ                         |



ʾ��ͼ������Ӧ�ó�������ͼ��ʾ��ͼ�͸�ʽ��һ�����������ͼ���Լ�example


###  4 Example
------------

</br>See more detailed usage in the&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b>[Chinese Documentation](https://github.com/hewm2008/PopGeneCNV/blob/main/PopGeneCNVʹ���ֲ�_manual_chinese.pdf)</b>
</br>See more detailed usage in the&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b>[English Documentation](https://github.com/hewm2008/PopGeneCNV/blob/main/PopGeneCNVʹ���ֲ�_manual_chinese.pdf)</b>
</br>See the example directory and  Manual.pdf for more detail.
</br>��������  Manual.pdf for more detail �����ʾ�����ݺͽű������ڽ���ĳЩ��ַ�ͷ�һЩ�̳�
</br></br> 
           ../bin/PopGeneCNV	-i	Bam.list	-g	~/01.Ref/MC.gff	-ooutPrefix -t 5
</br>  Ŀ¼  example/�����������������ͽű��÷���


* Example 1)ǧ��VCF�ز���SNP������
</br> bam���ݴ�gffֻ�г�ǰ���һЩ�У����鿴��ʽ�� �����������棬��ʽ����������ʾ

������
out.gene.stat.gz 	������Ʒ�����������Ⱥ͸��Ƕȵ�ͳ����Ϣ
```
#GeneID SampleName      Length  CoveredSite     TotalDepth      Coverage%       MeanDepth
MCh00g14009     LB      1500    1500    1583610 100.00  1055.74
MCh00g14007     LB      261     161     177751  61.69   681.04
MCh00g14005     LB      246     146     396     59.35   1.61
MCh00g14004     LB      534     135     486     25.28   0.91
MCh00g14002     LB      222     222     141666  100.00  638.14

```

������
out.GeneDepth.mat.gz	������Ʒ�������ȣ������������(readcount��

```
#GeneID LB      LG      LH
MCh00g04914     0.00    0.00    0.00
MCh00g04915     0.00    0.00    0.00
MCh00g04916     0.00    0.00    0.00
MCh00g04917     1066.08 2306.49 1100.18
MCh00g04918     1050.94 2157.08 1061.07
MCh00g04919     439.54  883.84  456.20
```



������
Out.GeneCNV.Raw.mat.gz	��׼����Ⱥ�����CNV�ľ���


```
#GeneID LogPi��CNV��̬��ָ�꣩   S2(����)      LB(��Ʒ��)      LG(��Ʒ��)      LH(��Ʒ��)
MCh01g00249     0.00    0.00    0.83    0.74    0.87
MCh01g00250     0.00    0.00    1.07    1.02    0.96
MCh01g00251     0.00    0.01    1.10    0.90    1.07
MCh01g00252     0.00    0.00    1.12    1.25    1.15
MCh01g00256     0.00    0.00    0.77    0.77    0.75
```



������
Out.GeneCNV.filter.mat.gz	��׼���󣬴���Ⱥ���̬CNV�ĺ�ѡ����
```
#GeneID LogPi��CNV��̬��ָ�꣩   S2(����)      LB(��Ʒ��)      LG(��Ʒ��)      LH(��Ʒ��)
MCh02g06373     2.38    0.13    0.00    0.76    0.00
MCh02g06524     0.17    0.11    0.37    1.02    0.28
MCh02g06525     0.13    0.13    0.42    1.18    0.41
MCh02g06552     0.89    0.27    0.07    1.17    0.07
MCh04g17936     0.30    0.26    0.38    1.40    0.24
MCh04g17937     0.22    0.25    0.37    1.44    0.40
```


������
Out.GeneCNV.vcf.mat.gz	��׼���󣬴���Ⱥ���̬CNV�ĺ�ѡ���� ��λ�÷�����
```
##fileformat=VCF4.1
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplation">
##INFO=<ID=END,Number=1,Type=String,Description="End of the Gene Site">
##INFO=<ID=LogPi,Number=1,Type=String,Description="the LogPi of This Gene">
##INFO=<ID=S2,Number=1,Type=String,Description="variance of this CNV">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	LB	LG	LH
MC02	12054869	MCh02g06373	N	<DEL>	.	.	END=12055081;LogPi=2.38;S2=0.13	GT	1/1	0/0	1/1
MC02	16603048	MCh02g06524	N	<DEL>	.	.	END=16603437;LogPi=0.17;S2=0.11	GT	1/1	0/0	1/1
MC02	16620545	MCh02g06525	N	<DEL>	.	.	END=16620860;LogPi=0.13;S2=0.13	GT	1/1	0/0	1/1
MC02	17900213	MCh02g06552	N	<DEL>	.	.	END=17900503;LogPi=0.89;S2=0.27	GT	1/1	0/0	1/1
MC04	12042994	MCh04g17936	N	<DEL>	.	.	END=12043973;LogPi=0.30;S2=0.26	GT	1/1	0/0	1/1
MC04	12076848	MCh04g17937	N	<DEL>	.	.	END=12078351;LogPi=0.22;S2=0.25	GT	1/1	0/0	1/1
MC04	12088387	MCh04g17938	N	<DEL>	.	.	END=12093579;LogPi=0.16;S2=0.21	GT	1/1	0/0	1/1
MC04	12100198	MCh04g17939	N	<DEL>	.	.	END=12101215;LogPi=0.21;S2=0.15	GT	1/1	0/0	1/1
MC04	13601191	MCh04g17989	N	<DEL>	.	.	END=13606780;LogPi=0.85;S2=0.12	GT	1/1	0/0	1/1
MC04	13619213	MCh04g17991	N	<DEL>	.	.	END=13619620;LogPi=2.34;S2=0.12	GT	1/1	0/0	1/1
```




###  5 Advantages

</br>�ٶȿ죬���ڴ�    fast speed, low memory
</br>��������    Simple and easy to use
</br>�ⰲװ   Free installation


###  6 Discussing
------------
- [:email:](https://github.com/hewm2008/PopGeneCNV) hewm2008@gmail.com / hewm2008@qq.com
- join the<b><i> QQ Group : 125293663</b></i>

######################swimming in the sky and flying in the sea #############################

