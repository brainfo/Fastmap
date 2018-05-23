#root dir /home/jianghong/FASTmap/raw_data 
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-d', action='store', dest='dirname',
        help='input folder name')
parser.add_argument('-f',action='append',dest='sample',
        help='input sample name')
parser.add_argument('r',action='store',dest='ref',
        help='input reference name with .fa')

paras = parser.parse_args()
current_path = os.getcwd()
source_path = '%s/input/%s' %(current_path, paras.dirname)
target_path = '%s/output/%s' %(current_path, paras.dirname)
ref = '%s/reference/%s/%s' %(current_path, paras.dirname,paras.ref)
samples = paras.sample

def prepareFolders():
    if not os.path.exists(sorce_path): 
        print '%s not exist!' %sorce_path
        return False
    if not os.path.exists(target_path):
        os.mkdir(target_path)
        print 'Create %s successfully!' %target_path
        print 'Results will write into %s' %target_path
    if not os.path.exists(target_path+'/reference'): os.mkdir(target_path+'/reference')
    return True

## =============================
## functions for create index
## Author: jianghong
## Create date: 2018-03-10
## =============================

def hasBWAIndex(ref):
    appendix = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    for app in appendix:
        if not os.path.isfile(ref + app):
            return False
    return True

def hasSamtoolsIndex(ref):
    appendix = ['.fai']
    for app in appendix:
        if not os.path.isfile(ref + app):
            return False
    return True

def hasPICARDIndex(ref):
    appendix = ['.dict']
    ref=ref[:-3]
    for app in appendix:
        if not os.path.isfile(ref + app):
            return False
    return True

def createIndex(ref):
    ## step1: build bwa index
    if not hasBWAIndex(ref):
        cmd = 'bwa index %s' %(ref)
        os.system(cmd)
    ## step2: build samtools index
    if not hasSamtoolsIndex(ref):
        cmd = 'samtools faidx %s' %(ref)
        os.system(cmd)    
    ## step3: build picard index
    if not hasPICARDIndex(ref):
        cmd = 'java -jar CreateSequenceDictionary.jar R=%s O=%s.dict' %(ref, ref)
        os.system(cmd)

## ==============Update: 2018-03-10===================
## ==============By: JarningGau ======================

## ================================
## Functions for variants calling
## ================================


def alignsample(sample,ref):
    outfile = target_path + sample + '.bam'
    if not os.path.isfile(outfile):
        read1 = source_path + sample + '_1.clean.fq.gz'
        read2 = source_path + sample + '_2.clean.fq.gz'
        cmd = 'bwa mem %s %s %s | samtools view -q 30 -S -b - > %s' %(ref, read1, read2, outfile)
        os.system(cmd)
    return outfile

def sortbam(sample,ref):
    outfile = target_path + sample + '.sorted.bam'
    if not os.path.isfile(outfile):
        infile = alignsample(sample,ref)
        cache_bam = target_path + sample + '.nnnn.bam'
        cmd = 'samtools sort -T %s -O bam -o %s %s' %(cache_bam, outfile, infile)
        os.system(cmd)
    return outfile

def markduplicates(sample,ref):
    outfile = target_path + sample + '.sorted.markdup.bam'
    if not os.path.isfile(outfile):
        infile = sortbam(sample,ref)
        print infile
        metrics = target_path + sample + '.markdup_metrics.txt'
        cmd = 'java -jar /home/public/tools/picard-tools-1.124/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=%s O=%s M=%s' %(infile,outfile,metrics)
        os.system(cmd)
    baifile = outfile + '.bai'
    if not os.path.isfile(baifile):
        infile = outfile
        cmd = 'samtools index %s' %(infile)
        os.system(cmd)
    return outfile


def callvariants(sample,ref):
    outfile = target_path + sample + '.vcf.gz'
    if not os.path.isfile(outfile):
        infile = markduplicates(sample,ref)
        cmd = 'samtools mpileup -ugf %s %s -C 50 -d 100| bcftools call -vmO z -o %s' %(ref, infile, outfile)
        os.system(cmd)
        cmd = 'vcftools --gz %s --remove-indels --recode --recode-INFO-all --out %s_SNPs_only' %(outfile, sample)

## ================================

if __name__ == '__main__':
    ## creat index
    if not prepareFolders(): return 0
    createIndex(ref)
    for sample in samples:
        callvariants(sample,ref)