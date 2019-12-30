#! python

import argparse
from Bio import AlignIO,SeqIO
import sys
import subprocess as subp
import os
import regex

parser=argparse.ArgumentParser(description='This program does alignment of protein or DNA sequences using muscle, trim the alignment using Gblock, and reconstruct ML phylogenetic trees using RAxML or FastTree. Do NOT run multiple sequence files in the same folder because they generate the same intermediate files.')
parser.add_argument('-i','--input',dest='inp',required=True,help='Input fasta file')
parser.add_argument('-o','--output',dest='out',required=False,default='ml.nwk',help='Output newick file of tree [Default: ml.nwk]')
parser.add_argument('-m','--maker',dest='maker',required=False,default='RAxML',help='Tree maker, FastTree or RAxML [Default: RAxML]')
parser.add_argument('-t','--type',dest='type',required=False,default='d',help='Data type, d for DNA, p for protein, c for codon [Default: d]')
parser.add_argument('-s','--model',dest='model',required=False,default='default',help='Substitution model, RAxML: GTRCAT, GTRCATI, GTRGAMMA or GTRGAMMI for nucleotide data, PROTGAMMALG or PROTGAMMAJTT for amino acid data; FastTree: \'-gtr -gamma\' for nucleotide data, \'-lg -gamma\' for amino acid data [Default: RAxML => GTRGAMMA for nucleotide, PROTGAMMALG for amino acid; FastTre => -gtr -gamma for nucleotide, -lg -gamma for amino acid]')
parser.add_argument('-n','--bootstrap',dest='boot',required=False,default='100',help='Number of bootstrap replicates, only for RAxML [Default: 100]')
parser.add_argument('-c','--cpu',dest='cpu',required=False,default='10',help='CUP cores to use [Default: 10]')

args=parser.parse_args()

home=os.path.dirname(os.path.realpath(__file__))+'/bin'

# muscel alignment
if not subp.call('set -ex; '+home+'/muscle -in '+args.inp+' -out '+args.inp+'.aln',shell=True)==0:
    print('Error: Alignment failed.')
    sys.exit(1)
else:
    print('Alignment done.')

# rename alignment
outfile=open(args.inp+'.aln.renamed','w')
num=0
di={}
for rec in SeqIO.parse(args.inp+'.aln','fasta'):
    num+=1
    di[str(num)]=rec.id
    rec.id=str(num)
    rec.description=''
    SeqIO.write(rec,outfile,'fasta')

outfile.close()
os.system('mv '+args.inp+'.aln.renamed '+args.inp+'.aln')
print('Rename the alignment done.')

# gblock trimming
subp.call('set -ex; '+home+'/gblocks '+args.inp+'.aln -t='+args.type+' -b4=3 -b5=h',shell=True)
subp.call('rm -rf '+args.inp+'.aln-gb.htm',shell=True)

if not os.path.isfile(args.inp+'.aln-gb'):
    print('Error: Gblocks trimming failed.')
    sys.exit(1)
else:
    print('Gblocks trimming done.')

# convert to phylip
if not subp.call('set -ex; python '+home+'/fasta2phylip.py -i '+args.inp+'.aln-gb -o '+args.inp+'.aln-gb.phy',shell=True)==0:
    print('Error: Converting to phylip failed.')
    sys.exit(1)
else:
    print('Conversion to Phylip done.')

# raxml
if args.maker=='RAxML':
    os.system('rm -rf RAxML_*')
    model={'d':'GTRGAMMA','p':'PROTGAMMALG','c':'GTRGAMMA'}

    if args.model=='default':
        args.model=model[args.type]

    if not subp.call('set -ex; '+home+'/raxml -f a -k -s '+args.inp+'.aln-gb.phy -m '+args.model+' -n '+args.out+' -x 12345 -p 12345 -N '+args.boot+' -T '+args.cpu,shell=True)==0:
        print('Error: RAxML running failed.')
        sys.exit(1)
    else:
        print('RAxml running done.')

    os.system('mv RAxML_bipartitionsBranchLabels.'+args.out+' '+args.out)
    os.system('rm -rf RAxML_*')

elif args.maker=='FastTree':
    model={'d':'-gtr -gamma','p':'-lg -gamma','c':'-gtr -gamma'}

    if args.model=='default':
        args.model=model[args.type]
    if args.type=='d' or args.type=='c':
        args.type='-nt'
    elif args.type=='p':
        args.type=''
    else:
        print('Error: incorrect data type')
        sys.exit(1)

    if not subp.call('set -ex; '+home+'/fasttree -out '+args.out+' '+args.model+' '+args.type+' '+args.inp+'.aln-gb.phy',shell=True)==0:
        print('Error: FastTree running failed.')
        sys.exit(1)
    else:
        print('FastTree running done.')

else:
    print('Error: incorrect maker type')
    sys.exit(1)


infile=open(args.out,'r')
tree=infile.read()
tree=regex.sub(r'([,\(])(\d+):',lambda x:x.group(1)+di[x.group(2)]+':',tree)
infile.close()
outfile=open(args.out,'w')
outfile.write(tree)
outfile.close()





