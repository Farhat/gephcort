#!/usr/bin/python

'''

####################################################################################
# This short script is a part of gephcort package v1.0
# 
# For any queries related to script, contact: Amol Kolte (amolkolte1989@gmail.com)
####################################################################################

Created on 29-Mar-2012

'''

####################################################################################
# Scanning Inputs
####################################################################################

import sys, getopt, copy, pickle
from multiprocessing import Pool
import time
import redis

# initializing redis server
REDIS = redis.Redis("localhost")


def main(argv):
    global seq, intree, seq_format, iterations, phen, out, ressurect_file, log_file, num_cores

    try:
        opts, args = getopt.getopt(argv, "hs:t:f:i:p:o:r:l:c:", ['seq=', 'intree=', 'seq_format=', 'iterations=', 'phen=', 'out=', 'ressurect_file=', 'log_file=', "cores="])
    except getopt.GetoptError:
        print "\n USAGE: python reanimate.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -r <ressurect_output_file> -o <output_file> \n \n \t --seq, -s : SNP sequence file \n \t --tree, -t : Newick tree \n \t --seq_format, -f : SNP sequence file format (phylip/fasta) \n \t --iter, -i : Phenotype shuffling iterations \n \t --phen, -p : Custom format phenotype file \n \t --out, -o : Output filename \n \t --ressurect_file, -r : Output file obtained from ressurect.R for a given seq and tree file \n \t --log_file, -l : Log file (Optional)\n"

    for opt, arg in opts:
        if opt=='-h':
            print '\n\n python reanimate.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -r <ressurect_output_file> -o <output_file> -l <log_file> (Optional)\n\n'
            sys.exit()        
        elif opt in ("-s", "--seq"):
            seq=arg
        elif opt in ("-t", "--intree"):
            intree=arg
        elif opt in ("-f", "--seq_format"):
            if arg=="phylip":
                seq_format='iphylip'
            elif arg=="fasta":
                seq_format='fasta'
            else:
                print "\n** Fatal Error : Sequence format not supported !!! Kindly consider converting it into 'fasta' or 'phylip' format **"
        elif opt in ("-i", "--iterations"):
            iterations=arg
        elif opt in ("-p", "--phen"):
            phen=arg
        elif opt in ("-o", "--out"):
            out=arg
        elif opt in ("-r", "--ressurect_file"):
            ressurect_file=arg
        elif opt in ("-l", "--log_file"):
            log_file=arg
        elif opt in ("-c", "--cores"):
        	num_cores = int(arg)
            

if __name__ == "__main__":
    if sys.argv[1:]:
        main(sys.argv[1:])
    else:
        print "\n USAGE: python reanimate.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -r <ressurect_output_file> -o <output_file> \n \n \t --seq, -s : SNP sequence file \n \t --tree, -t : Newick tree \n \t --seq_format, -f : SNP sequence file format (phylip/fasta) \n \t --iter, -i : Phenotype shuffling iterations \n \t --phen, -p : Custom format phenotype file \n \t --out, -o : Output filename \n \t --ressurect_file, -r : Output file obtained from ressurect.R for a given seq and tree file \n \t --log_file, -l : Log file (Optional) \n"
        sys.exit()


###############################################################################
# Loading required modules
###############################################################################
time_elapsed = {}
start_time = time.time()
ssss = time.time()
print "Initializing .. Loading modules"
if num_cores == None:
	num_cores = 1

import time, scipy, math, numpy as np # system commands
from ete2 import PhyloTree    # To create a phylogenetic tree
from random import shuffle      # Shuffle phenotype
from scipy import stats         # To calculate p-value from z-score
from rpy2.robjects.packages import importr
from sets import Set


# Opening log file
try:
    log=open(log_file, "w")
except:
    log=open("gephcort_run.log", "w")
# Logging start time
log.write("Start time: "+str(time.localtime()[0])+"-"+str(time.localtime()[2])+"-"+str(time.localtime()[1])+"\t"+str(time.localtime()[3])+":"+str(time.localtime()[4])+":"+str(time.localtime()[5])+"\n")

ape_objects={"delta.plot":"delta_plot", "dist.dna":"dist_dna", "dist.nodes":"dist_nodes", "node.depth":"node_depth", 
"node.depth.edgelength":"node_depth_edgelength","node.height":"node_height", "node.height.clado":"node_height_clado", 
"prop.part":"prop_part"}

ape=importr("ape", robject_translations = ape_objects)    # Required for phangorn
ph=importr("phangorn")    # Phylogenetic operations in R

print "All modules imported successfully"

t = PhyloTree(intree, alignment=seq, alg_format=seq_format)     # Main tree containing entire sequence
dtp = PhyloTree(intree)                 # Dummy tree for phenotype shuffling

print "Tree file read successfully"

phenfile=open(phen, "r")    # Phenotype file
phenlist=[]
for line in phenfile.readlines():
    phenlist.append([line.split("\t")[0].strip(), line.split("\t")[1].strip()])
phenfile.close()

phenotype={}    # Dictionary containing species names and their phenotype values

# Phenotype file should have two columns separated by tab containing taxa name
# in the first column and a numerical phenotype value in the second
#
# dog1  3.2
# dog2  4.4
# cat2  4.5
# .
# .


if len(phenlist[1])==3:
    for i in range(len(phenlist)):
        try:
            phenotype[phenlist[i][0]] = [float(phenlist[i][1])-float(phenlist[i][2]), float(phenlist[i][1])+float(phenlist[i][2])]
        except IOError:
            print "Error while readling phenotype file"
if len(phenlist[1])==2:
    for i in range(len(phenlist)):
        try:
            phenotype[phenlist[i][0]] = [float(phenlist[i][1]), float(phenlist[i][1])]
        except IOError:
            print "Error while readling phenotype file"
else:
    print "Phenotype file format not supported, please go through the instructions"

missing = 18    # phangorn value for missing genotype
shphenotype=[]
shphenotype.append(phenotype.values())

#########################################################################################
# Reading complete ancestral sequence data generated through R (In the form of pseudo-patterns)
#########################################################################################
begin_tm = time.time()
print "Reading complete ancestral sequence data generated through R " + str(time.time() - start_time)
start_time = time.time()

rtree=PhyloTree(intree)    # Tree for "R" generated patterns
tree=ape.read_tree(intree)

rlist=[]
ropf = open(ressurect_file, "r")    # rdata.dat is a rgp.R generated output file
for tab in ropf.readlines():
    tab=tab.rstrip()
    rlist.append(tab.split(" "))
ropf.close()

ori=np.array(rlist)

for node in rtree.traverse("postorder"):    # Patterns are being linked to their corresponding nodes
    if node.is_leaf():
        node.add_features(data=[None for i in range(len(rlist[0])-1)])  # Its rlist[0]-1, because nucleotides begins with name of species
        for i in range(len(ori[:,0])):
            if '"'+node.name+'"' == ori[:,0][i] :
                node.add_features(rtoken=i+1)
    else :
        node.add_features(data=[None for i in range(len(rlist[0])-1)])
        node.add_features(rtoken=None)

for node in rtree.traverse("postorder"):
    if node.is_leaf():
        node.data=map(lambda x: x, rlist[node.rtoken-1][1:])    # Its the sequence after name
        node.up.rtoken=int(ph.Ancestors(tree, node.rtoken, "parent")[0])
    else:
        try:
            node.data=map(lambda x: x, rlist[node.rtoken-1][1:])
            node.name=node.rtoken
            node.up.rtoken=int(ph.Ancestors(tree, node.rtoken, "parent")[0])
        except :
            print "Root node is encountered, ancestral node mapping complete"

###############################################################################################
# Phenotype manipulations
###############################################################################################
print "Starting Phenotype manipulations " + str(time.time() - start_time)
start_time = time.time()

                              
def phrange(node):
    ''' Recursive function '''

    left, right = node.children
    if left.phenrange is None:
        phrange(left)
    if right.phenrange is None:
        phrange(right)
    
    temp = []
    nwrange = [None, None]
    temp.append(left.phenrange)
    temp.append(right.phenrange)
    temp.sort()

    if temp[1][0] < temp[0][1]:
        nwrange[0] = temp[1][0]
        nwrange[1] = min(temp[0][1], temp[1][1])
    elif temp[1][0] > temp[0][1]:
        nwrange[0] = temp[0][1]
        nwrange[1] = temp[1][0]
    else:
        nwrange[0] = temp[1][0]
        nwrange[1] = temp[0][1]
    node.phenrange = nwrange
    node.phenvalue = np.average(node.phenrange)
    left.phenvalue = np.average(left.phenrange)
    right.phenvalue = np.average(right.phenrange)



counter = 0

for node in dtp.traverse("preorder"):
    node.add_features(counter=counter)
    counter += 1
    if node.is_leaf():
        node.add_features(phenrange=phenotype[node.name], phenvalue=None)
    else:
        node.add_features(phenrange=None, phenvalue=None) 


pickle.dump(dtp, open("dtp.p","wb"), pickle.HIGHEST_PROTOCOL)



# List to contain phenvalues for all the nodes after several shuffling attempts.

sflphen = []        # List to contain phenvalues for all the nodes after several shuffling attempts.

phrange(dtp)

for p in xrange(int(iterations)):
    shuffle(shphenotype[0])
    sflphen.append([])

    index = 0   
    for node in dtp.traverse("preorder"):     # serves two purpose, 1) Add existing values to sflphen 2) reset tree phenrange attribites
        sflphen[p].append(node.phenvalue)
        if node.is_leaf():
            node.phenrange = shphenotype[0][index]
            index += 1 
        else:
            node.phenrange = None
            
    phrange(dtp)

##########################################################################################
# Generating patterns from R output data
##########################################################################################
print "Generating patterns from R output data"

single=True
for node in rtree.traverse("postorder"):
    if single==True:
        if node.is_leaf():
            ref=node
            single=False

dt =  rtree

pattern = {}    # Contains all the possible patterns of SNPs across the tree, generated by rgp.R

sfldict = {}    # Details of branches where nucleotide is changing and the corresponding p-value

bchanges = {}    # Number of braches the nucleotide is changing and the corresponding pattern

#mafreq={}    # Stores minor allele frequency for individal SNP

counter = 0
for node in dt.traverse('preorder'):
    node.add_features(counter=counter)
    counter += 1
############################################################################################################################################################################################################
# adding the values to be processed in the Queue.
ptns = set([])
total_queue_size = 0
for var in xrange(len(ref.data)):

    currentdata = []
    currentdict = {"18" : 18}

    token = 0
    for node in dt.traverse("postorder"):
       
        if node.is_leaf():
            try:
                currentdata.append(currentdict[node.data[var]])        # Pattern recognition
            except:
                currentdata.append(token)
                currentdict[node.data[var]] = token
                token += 1                
        else:
            pass   
    tdata = tuple(currentdata)

    total_queue_size += 1
###########################################################################################################################################################################################################

# prefix dict for lookup in redis.

prefix_dict = {}
prefix_dict['pattern'] = 'pat'
prefix_dict['sfldict'] = 'sfl'
prefix_dict['bchanges'] = 'bch'

tests = 0

def compute(input_data):
    tdata = str(input_data[0])
    temp = input_data[1]
    key = prefix_dict['pattern'] + tdata
    if REDIS.get(key) == None:
        temp.sort()
        tpp = map(tuple, temp)
        tp = str(tuple(tpp))
        dtemp = []
        sfl_key = prefix_dict['sfldict'] + tp
        sfl_value = REDIS.get(sfl_key)
        if sfl_value == None:
            for index in xrange(len(sflphen)):
                tsum = 0.0
                if not len(temp) == 0:
                    for n in xrange(len(temp)):
                        tsum += abs(sflphen[index][temp[n][0]]-sflphen[index][temp[n][1]])
                    dtemp.append(math.pow((tsum/len(temp)),1))
                else:
                    dtemp.append('NA')  # if all alleles are identical

            if 'NA' in dtemp:
                pvalue = -1
            else:
                pvalue = scipy.stats.norm.sf(abs(dtemp[0] - np.average(dtemp[1:])) / np.std(dtemp[1:])) * 2
                #pvalue=scipy.stats.ttest_1samp(dtemp[1:], dtemp[0])[1]
            
            REDIS.set(prefix_dict['sfldict'] + tp, pvalue)
        else:
            pvalue = float(sfl_value)
        REDIS.set(prefix_dict['pattern'] + tdata, pvalue)
        REDIS.set(prefix_dict['bchanges'] + tdata, len(tp))

temp_ray = []

for x in range(total_queue_size):
    temp = []
    currentdata = []
    currentdict = {"18" : 18}

    token = 0
    for node in dt.traverse("postorder"):
       
        if node.is_leaf():
            try:
                currentdata.append(currentdict[node.data[x]])        # Pattern recognition
            except:
                currentdata.append(token)
                currentdict[node.data[x]] = token
                token += 1                
        else:
            pass   
    tdata = tuple(currentdata)
    for node in dt.traverse('preorder'):       
            if not node.is_leaf():             
                left, right = node.children   
                if left.data[x] != node.data[x]:         # Nucleotide substitution
                    if left.data[x] != 18 :
                        temp.append([node.counter, left.counter])
                if right.data[x] != node.data[x]:         # Nucleotide substitution
                    if right.data[x] != 18 :
                        temp.append([node.counter, right.counter])
    
    temp_ray.append((tdata, temp))

p = Pool(num_cores)

p.map(compute, temp_ray)

p.close()
p.join() 

###############################################################################################
# Reading orignal sequence, evaluating by breaking into patterns
###############################################################################################
print "Reading orignal sequence, evaluating by breaking into patterns " + str(time.time() - start_time)
start_time = time.time()


single=True

for node in t.traverse("postorder"):
    if single == True:
        if node.is_leaf():
            oriref = node
            single = False

plist=[None for var in xrange(len(oriref.sequence))]
blist=[None for var in xrange(len(oriref.sequence))]
#maflist=[None for var in xrange(len(oriref.sequence))]

for var in xrange(len(oriref.sequence)):
    currentdata = []
    currentdict = {"-" : 18}
    token = 0
    for node in t.traverse("postorder"):
        if node.is_leaf():
            try:
                currentdata.append(currentdict[node.sequence[var]])        #Pattern recognition
            except:
                currentdata.append(token)
                currentdict[node.sequence[var]]=token
                token+=1                
        else:
            pass
              
    tdata = tuple(currentdata)
    
    try:
        
        plist[var] = float(REDIS.get(prefix_dict['pattern'] + str(tdata)))
        blist[var] = int(REDIS.get(prefix_dict['bchanges'] + str(tdata)))
#        maflist[var]=mafreq[tdata]
    except:
        
        print "Patterns missing"

#############################################################################################
# Performing multiple correction
#############################################################################################
uniq_plist=Set(plist)

try:
    uniq_plist.remove(-1)
except:
    pass

uniq_plist=list(uniq_plist)

sorted_uniq_plist=sorted(uniq_plist)

corrected={}    # corrected[raw_pvalue]=[fdr_pvalue, bonferroni_pvalue]

for entry in sorted_uniq_plist:
    corrected[entry]=[min(1, entry*len(sorted_uniq_plist)/(len(sorted_uniq_plist)-sorted_uniq_plist.index(entry))), min(1, entry*len(sorted_uniq_plist))]

#############################################################################################
# Writing output
#############################################################################################

try:
    opf=open(out, "w")
except:
    opf=open("/tmp/gpresult.dat", "w")

opf.write('#SNP_index'+'\t'+'p-value'+'\t'+'p.adjusted_FDR'+'\t'+'p.adjusted_Bonferroni'+'\n')
for i in range(len(plist)):
    if plist[i] in corrected:
        opf.write(str(i)+'\t'+str(plist[i])+'\t'+str(corrected[plist[i]][0])+'\t'+str(corrected[plist[i]][1])+'\n')
    else:
        opf.write(str(i)+'\t'+'NA'+'\t'+'NA'+'\t'+'NA'+'\n')
opf.close()

#############################################################################################
# Writing Summary Log
#############################################################################################

log.write("End time: "+str(time.localtime()[0])+"-"+str(time.localtime()[2])+"-"+str(time.localtime()[1])+"\t"+str(time.localtime()[3])+":"+str(time.localtime()[4])+":"+str(time.localtime()[5])+"\n\n")
log.write("Genotype file\t"+seq+"\n")
log.write("Newick file\t"+intree+"\n")
log.write("Phenotype file\t"+phen+"\n")
log.write("Iterations for permutation test\t"+str(iterations)+"\n")
log.write("Number of SNPs\t"+str(len(oriref.sequence))+"\n")
log.write("Total number of observed patters\t"+ str(len(pattern))+"\n")
log.write("Number of statistical tests performed\t"+str(tests))
log.close()
print "Operation complete", time.time() - ssss

#################################################################################################

