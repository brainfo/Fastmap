#coding=utf-8
# ==== FASTmap is using for identify candidate mutations ====
# ==== based on linkage analysis and site filtering in   ====
# ==== genetic screening of C.elegans                    ====
# #### Date: 2018-06-01 ####
# #### Version: 1.0 ####

import os
import argparse
import pandas as pd 
import numpy as np 
import re 


## Parameters
parser = argparse.ArgumentParser()
parser.add_argument('-a', action='store_const', dest='mode',
        const='annotation',
        help='annotation mode')
parser.add_argument('-l', action='store_const', dest='mode',
        const='linkage',
        help='linkage_analysis mode')
parser.add_argument('--window', action='store', dest='window_size', default=0.4, 
        help='window size for smooth', type=float)
parser.add_argument('--slide', action='store', dest='step_size', default=0.1, 
        help='step size for smooth', type=float)
parser.add_argument('-d', action='store', dest='dirname',
        help='input folder name')
parser.add_argument('--chr', action='store', dest='chrom',
        help='select chromsome')
parser.add_argument('-b', action='store', dest='begin',
        help='select region begin', type=float)
parser.add_argument('-e', action='store', dest='end',
        help='select region end', type=float)
parser.add_argument('-c', action='append',dest='depth',
        help='normal depth range', type=int)
parser.add_argument('--version', action='version', version='FASTmap 1.0')


paras = parser.parse_args()
current_path = os.getcwd()
sorce_path = '%s/output/%s' %(current_path, paras.dirname)
target_path = '%s/R_result/%s' %(current_path, paras.dirname)
annovar_path = '%s/tools/annovar' 

def prepareFolders():
    if not os.path.exists(sorce_path): 
        print '%s not exist!' %sorce_path
        return False
    if not os.path.exists(target_path):
        os.mkdir(target_path)
        print 'Create %s successfully!' %target_path
        print 'Results will write into %s' %target_path
    if not os.path.exists(target_path+'/picture'): os.mkdir(target_path+'/picture')
    return True


## This function is using for put mut.vcf file into hash.
## Only hom records is put into this hash
def depth(line,low=10,high=100):
    #depth_range = paras.depth
    line_list = line.strip('\n').strip('\r').split('\t')
    info = line_list[7].split(';')
    DP4 = info[-2]
    pattern = re.compile(r'\d+')
    par = pattern.findall(DP4)
    depth = int(par[1]) + int(par[2]) + int(par[3]) + int(par[4])
    if depth > low & depth< high: return True
    else: return False

def putIntoHash(wt_file):
    h = {}
    with open(wt_file, 'r') as fi:
        for line in fi:
            if '#' in line: continue
            if '1/1' not in line: continue
            line_list = line.strip('\n').strip('\r').split('\t')
            mut = '%s_%s_%s_%s' %(line_list[0], line_list[1], line_list[3], line_list[4])
            # mut = I_11241_ACC_A
            h[mut] = 1
    return h


# Strip hom both exist in wt and mut
def tripWT(wt_file, mut_file, markers):
    h_WT = putIntoHash(wt_file)
    with open(markers, 'w') as fo:
        with open(mut_file, 'r') as fi:
            for line in fi:
                if '#' in line: continue
                #if not depth(line): continue
                line_list = line.strip('\n').strip('\r').split('\t')
                mut = '%s_%s_%s_%s' %(line_list[0], line_list[1], line_list[3], line_list[4])
                if mut not in h_WT: 
                    fo.write(line)


def R_ValueFunc(ref, alt):
    ref = float(ref)
    alt = float(alt)
    if ref == 0: return alt-1
    if alt == 0: return ref-1
    if ref == alt: return 0.
    if ref < alt: return alt/ref
    if ref > alt: return ref/alt
    #return (alt+1)/(alt+ref+2)


def smooth(data, start, end):
    '''
    data -- pandas.DataFrame
    start -- int
    end -- int
    return (pos_average, value_average)
    '''
    data_screen = data.loc[(data['Pos_Ref(Mb)'] > start) & (data['Pos_Ref(Mb)'] < end)]
    if data_screen.shape[0] < 1: return 0,0
    pos = np.average(data_screen['Pos_Ref(Mb)'])
    value = np.average(data_screen.R_value)
    return pos,value


def prePlot(data, chrom, win_size = 0.3, slide = 0.1):
    with open('%schrom_%s_for_plot.txt' %(target_path, chrom), 'w') as fo:
        data_chrom = data.loc[data.Chrom == chrom]
        y_max=data_chrom.R_value.max()
        slide_times = int(data['Pos_Ref(Mb)'].max() / slide)
        if data_chrom.depth>22 and data_chrom.depth<96 and data_chrom.QD>2:
            for i in xrange(slide_times + 1):
                start = i * slide
                end = start + win_size
                pos,value = smooth(data_chrom, start, end)
                if (pos,value) == (0,0): continue
                fo.write('%.4lf\t%.4lf\n' %(pos, value))
    return y_max


def extractMarkers():
    wt_file = '%swt.vcf' %source_path
    mut_file = '%smut.vcf' %source_path
    markers = '%smarkers.vcf' %target_path
    tripWT(wt_file, mut_file, markers)


def calculateR_value():
    # Input and output
    markers = '%smarkers.vcf' %target_path
    result_file = '%smut_R.txt' %target_path
    # ====
    with open(result_file, 'w') as fo:
        fo.write("Chrom\tPos_Abo(nt)\tPos_Ref(Mb)\tRO,AO\tR_value\n")
        with open(markers, 'r') as fi:
            for line in fi:
                line_list = line.strip('\n').strip('\r').split('\t')
                info = line_list[7].split(';')
                DP4 = info[-2]
                pattern = re.compile(r'\d+')
                par = pattern.findall(DP4)
                ref = int(par[1]) + int(par[2])
                alt = int(par[3]) + int(par[4])
                depth = ref +alt
                QD = line_list[5]/depth
                r_value = R_ValueFunc(ref, alt)
                Pos_Abo = line_list[1]
                Pos_Ref = int(Pos_Abo)/1e6
                fo.write("%s\t%s\t%.3lf\t%d,%d\t%.3lf,%d\t%.3lf\n" %(line_list[0], Pos_Abo, Pos_Ref, ref, alt, r_value,depth,QD))


def plotMain(win_size, slide):
    fileIn = '%smut_R.txt' %target_path
    data = pd.read_csv(fileIn, sep='\t')
    for chrom in ['I', 'II', 'III', 'IV', 'V', 'X']:
        y_max=prePlot(data, chrom, win_size, slide)
    for chrom in ['I', 'II', 'III', 'IV', 'V', 'X']:
        file = 'chrom_%s_for_plot.txt' %chrom
        img = 'chrom_%s.pdf' %chrom
        cmd = "Rscript /home/jianghong/FASTmap/script/plot.R %s %s %s %s" %(file, img, y_max, target_path)
        os.system(cmd)
    #with open('%spicture/parameters.txt' %target_path, 'w') as fo:
        #fo.write('window_size=%s\nstep_size=%s\n' %(win_size, slide))


def selectRegion(chrom, begin, end):
    if not os.path.exists(target_path+'/annotation'): os.mkdir(target_path+'/annotation')
    begin = int(begin * 1e6)
    end = int(end * 1e6)
    with open('%s/annotation/candidates.txt' %target_path, 'w') as fo:
        with open('%s/markers.vcf' %target_path, 'r') as fi:
            #base = ['A','C']
            for line in fi:
                line_list = line.strip('\n').strip('\r').split('\t')
                if line_list[0] != chrom: continue
                #if len(line_list[3]) != len(line_list[4]): continue # small indels are not EMS-like mutation
                if '1/1' not in line: continue  # only hom is selected
                if not depth(line,10,100): continue
                if int(line_list[1]) < begin: continue
                if int(line_list[1]) > end: continue
                fo.write('%s\t%s\t%d\t%s\t%s\t%s\n' 
                    %(line_list[0], line_list[1], int(line_list[1])+len(line_list[3])-1, line_list[3], line_list[4], line_list[5]))


def annotation():
    annovar_path = '%s/tools/annovar' %current_path
    fileIn = '%s/annotation/candidates.txt' %target_path
    cmd = '%s/table_annovar.pl %s %s/celdb/ --outfile %s/annotation/candidates --buildver cel --protocol refGene --operation g' %(
        annovar_path, fileIn, annovar_path, target_path)
    os.system(cmd)

def linkage_main():
    if not prepareFolders(): return 0
    extractMarkers()
    caculateR_value()
    plotMain(paras.window_size, paras.step_size)


def annotation_main():
    #selectRegion(paras.chrom, paras.begin, paras.end)
    selectRegion('III', 4, 10)
    annotation()


if __name__ == '__main__':
    #if paras.mode == 'annotation':
        annotation_main()
    #elif paras.mode == 'linkage':
        #linkage_main()