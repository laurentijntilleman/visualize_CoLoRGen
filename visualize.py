#!/usr/bin/env python3

# packages

import gzip
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging

class VisualizeHaplotypes:

    def __init__(self,start,stop,chromosome,variant_ref,gene,qual=True,variant_type=False):
        self.start = start
        self.stop = stop
        self.chromosome = chromosome
        self.gene = gene
        self.qual = qual
        self.variant_type = variant_type
        [self.haplo1,self.haplo2] = self.get_variants(variant_ref)

    def get_variants(self,variant_ref):
        haplo1 = []
        haplo2 = []
        nucleotides = ['A','T','G','C']
        if variant_ref.endswith('.gz'):
            vcf_reader = gzip.open(variant_ref)
        else:
            vcf_reader = open(variant_ref)
        for line in vcf_reader.readlines():
            if type(line) == str:
                l2 = line
            else:
                l2 = line.decode("utf-8")
            if not l2.startswith('#'):
                variant = l2.split()
                if self.qual:
                    if variant[6] == 'PASS':
                        if variant[0] == self.chromosome:
                            if int(variant[1]) >= int(self.start) and int(variant[1]) <= int(self.stop):
                                if variant[3] in nucleotides and variant[4] in nucleotides:
                                    if self.variant_type != 'INDEL':
                                        if '/' in variant[-1].split(':')[0]:
                                            haplo = variant[-1].split(':')[0].split('/')
                                        if '|' in variant[-1].split(':')[0]:
                                            haplo = variant[-1].split(':')[0].split('|')
                                        if int(haplo[0]) > 0:
                                            haplo1.append(int(variant[1]))
                                        if int(haplo[1]) > 0:
                                            haplo2.append(int(variant[1]))
                                elif self.variant_type != 'SNP':
                                    if '/' in variant[-1].split(':')[0]:
                                        haplo = variant[-1].split(':')[0].split('/')
                                    if '|' in variant[-1].split(':')[0]:
                                        haplo = variant[-1].split(':')[0].split('|')
                                    if int(haplo[0]) > 0:
                                        haplo1.append(int(variant[1]))
                                    if int(haplo[1]) > 0:
                                        haplo2.append(int(variant[1]))
                else:
                    if variant[0] == self.chromosome:
                        if int(variant[1]) >= int(self.start) and int(variant[1]) <= int(self.stop):
                            if variant[3] in nucleotides and variant[4] in nucleotides:
                                if self.variant_type != 'INDEL':
                                    if '/' in variant[-1].split(':')[0]:
                                        haplo = variant[-1].split(':')[0].split('/')
                                    if '|' in variant[-1].split(':')[0]:
                                        haplo = variant[-1].split(':')[0].split('|')
                                    if int(haplo[0]) > 0:
                                        haplo1.append(int(variant[1]))
                                    if int(haplo[1]) > 0:
                                        haplo2.append(int(variant[1]))
                            elif self.variant_type != 'SNP':
                                if '/' in variant[-1].split(':')[0]:
                                    haplo = variant[-1].split(':')[0].split('/')
                                if '|' in variant[-1].split(':')[0]:
                                    haplo = variant[-1].split(':')[0].split('|')
                                if int(haplo[0]) > 0:
                                    haplo1.append(int(variant[1]))
                                if int(haplo[1]) > 0:
                                    haplo2.append(int(variant[1]))
        vcf_reader.close()
        return [haplo1,haplo2]

    def check_variants(self,ref_haplo,found_haplo):
        [haplo1,haplo2] = ref_haplo
        [haplo1_1,haplo2_1] = found_haplo
        haplo1_1_F = [x for x in haplo1_1 if x not in list(set(haplo1) | set(haplo2))]
        haplo2_1_F = [x for x in haplo2_1 if x not in list(set(haplo1) | set(haplo2))]
        haplo1_1_P1 = [x for x in haplo1_1 if x in list(set(haplo1) | set(haplo2))]
        haplo2_1_P1 = [x for x in haplo2_1 if x in list(set(haplo1) | set(haplo2))]
        a11 = len(set(haplo1_1_P1) & set(haplo1))
        a12 = len(set(haplo1_1_P1) & set(haplo2))
        a21 = len(set(haplo2_1_P1) & set(haplo1))
        a22 = len(set(haplo2_1_P1) & set(haplo2))
        a_max = max([a11,a12,a21,a22])
        if a_max == a11 or a_max == a22 :
            haplo1_1_T = list(set(haplo1_1_P1) & set(haplo1))
            haplo2_1_T = list(set(haplo2_1_P1) & set(haplo2))
            haplo1_1_P = list(set(haplo1_1_P1) - set(haplo1))
            haplo2_1_P = list(set(haplo2_1_P1) - set(haplo2))
            haplo1_1_F1 = haplo1_1_F
            haplo2_1_F1 = haplo2_1_F
        else:
            haplo1_1_T = list(set(haplo2_1_P1) & set(haplo1))
            haplo2_1_T = list(set(haplo1_1_P1) & set(haplo2))
            haplo1_1_P = list(set(haplo2_1_P1) - set(haplo1))
            haplo2_1_P = list(set(haplo1_1_P1) - set(haplo2))
            haplo1_1_F1 = haplo2_1_F
            haplo2_1_F1 = haplo1_1_F
        return [haplo1_1_F1,haplo2_1_F1,haplo1_1_T,haplo2_1_T,haplo1_1_P,haplo2_1_P]

    def plot_base(self,figsize=(10,10)):
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        for spine in ["left", "top", "right"]:
            ax.spines[spine].set_visible(False)
        ax.set_xlim(int(self.start),int(self.stop))
        ax.scatter(self.haplo1,np.zeros(len(self.haplo1)),marker=2,s=200,color='k')
        ax.scatter(self.haplo2,np.zeros(len(self.haplo2)),marker=3,s=200,color='k')
        plt.axhline(y=0, color='k')
        plt.xlabel(self.chromosome)
        return [plt,fig,ax]

    def plot_variant(self,variant_vcf,plot1,t=0):
        plt, fig, ax = plot1
        ref_haplo = [self.haplo1,self.haplo2]
        found_haplo = self.get_variants(variant_vcf)
        phased_haplo = self.check_variants(ref_haplo,found_haplo)
        [haplo1_1_F,haplo2_1_F,haplo1_1_T,haplo2_1_T,haplo1_1_P,haplo2_1_P] = phased_haplo
        ax.scatter(haplo1_1_F,np.ones(len(haplo1_1_F))*-t,marker=2,s=200,color='red')
        ax.scatter(haplo2_1_F,np.ones(len(haplo2_1_F))*-t,marker=3,s=200,color='red')
        ax.scatter(haplo1_1_P,np.ones(len(haplo1_1_P))*-t,marker=2,s=200,color='orange')
        ax.scatter(haplo2_1_P,np.ones(len(haplo2_1_P))*-t,marker=3,s=200,color='orange')
        ax.scatter(haplo1_1_T,np.ones(len(haplo1_1_T))*-t,marker=2,s=200,color='green')
        ax.scatter(haplo2_1_T,np.ones(len(haplo2_1_T))*-t,marker=3,s=200,color='green')
        plt.axhline(y=-t, color='k')

        return [plt, fig, ax]


    def plot(self,variant_vcf,labelsize = 20):
        figsize = (20,2)
        plot1 = self.plot_base(figsize)
        t=1
        plot1 = self.plot_variant(variant_vcf,plot1,t*2)
        t+=1
        [plt, fig, ax] = plot1
        ax.tick_params(axis = 'both', which = 'major', labelsize = labelsize)
        ax.xaxis.label.set_size(labelsize)
        ax.ticklabel_format(useOffset=False, style='plain')
        xticks_start = self.start
        xticks_end = self.stop
        xticks_2 = int(xticks_start + (xticks_end - xticks_start)/3)
        xticks_3 = int(xticks_start + (xticks_end - xticks_start)*2/3)
        ax.set_xticks([xticks_start,xticks_2,xticks_3,xticks_end])
        plt.yticks([])
        ax.set_ylim([-t*2+0.5,1])
        logging.debug('fig')
        logging.debug(fig)
        self.fig = fig
        return fig

    def savefig(self,folder,name='dot_plot'):
        self.fig.savefig('{}{}.png'.format(folder,name))
        self.fig.savefig('{}{}.pdf'.format(folder,name))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='get second alignments')
    parser.add_argument('--log_level', metavar='log_level', dest='log_level',
                        type=str, help='log level')
    parser.add_argument('--log_file', metavar='log_file', dest='log_file',
                        type=str, help='log file')
    parser.add_argument('--vcf', metavar='vcf', dest='variant_path',
                        type=str, help='path to the called variants')
    parser.add_argument('-r','--ref', metavar='ref', dest='variant_ref',
                        type=str, help='path to the reference variants')
    parser.add_argument('--start', metavar='int', dest='start',
                        type=int, help='start target site')
    parser.add_argument('--stop', metavar='int', dest='stop',
                        type=int, help='stop target site')
    parser.add_argument('--name', metavar='name', dest='name',
                        type=str, help='locus name')
    parser.add_argument('-c','--chromosoom', metavar='chr', dest='chromosome',
                        type=str, help='chromosoom of the target region')
    parser.add_argument('-f','--file', metavar='file', dest='file',
                        type=str, help='path to save figure')
    args = parser.parse_args()

    # logging
    log_numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log_level)
    logging.basicConfig(level=log_numeric_level, filename=args.log_file,
                        format='%(asctime)s %(message)s')

    logging.info('start snps')

    vis_snp = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.name,variant_type='SNP')
    fig_snp = vis_snp.plot(args.variant_path)
    vis_snp.savefig(args.file,f'{args.name}_SNP')

    logging.info('start indels')

    vis_indel = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.name,variant_type='INDEL')
    fig_indel = vis_indel.plot(args.variant_path)
    vis_indel.savefig(args.file,f'{args.name}_INDEL')

    logging.info('start all')

    vis = VisualizeHaplotypes(args.start,args.stop,args.chromosome,args.variant_ref,args.name,variant_type=None)
    fig = vis.plot(args.variant_path)
    vis.savefig(args.file,args.name)
