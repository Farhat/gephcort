
'''
Created on 29-Mar-2012

@author: amol
'''
#print "Initializing .. Loading parameters"

import sys, getopt, scipy, math, numpy as np # system commands
from ete2 import PhyloTree    # To create a phylogenetic tree
from random import shuffle      # Shuffle phenotype
from scipy import stats         # To calculate p-value from z-score
from rpy2.robjects.packages import importr

###############################################################################
# Scanning Inputs
###############################################################################

def main(argv):
	global seq, intree, seq_format, iterations, phen, out

	try:
		opts, args = getopt.getopt(argv, "hs:t:f:i:p:o:", ['seq=', 'intree=', 'seq_format=', 'iterations=', 'phen=', 'out='])
	except getopt.GetoptError:
		print "\n USAGE: phenotype_reconstruction_and_correlation.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -o <output_file> \n \n \t --seq, -s : SNP sequence file [seq] \n \t --tree, -t : Newick tree [tree] \n \t --seq_format, -f : SNP sequence file format (phylip/fasta) [format] \n \t --iter, -i : Phenotype shuffling iterations [iter] \n \t --phen, -p : Custom format phenotype file [phen] \n \t --out, -o : Outputfile redirective (not recommanded) [out]"

	for opt, arg in opts:
		if opt=='-h':
			print 'phenotype_reconstruction_and_correlation.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -o <output_file>'
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
			

if __name__ == "__main__":
	if sys.argv[1:]:
		main(sys.argv[1:])
	else:
		print "\n USAGE: phenotype_reconstruction_and_correlation.py -s <seq_file> -t <tree_file> -f <format(fasta/phylip)> -i <phen_iterations> -p <phen_file> -o <output_file> \n \n \t --seq, -s : SNP sequence file [seq] \n \t --tree, -t : Newick tree [tree] \n \t --seq_format, -f : SNP sequence file format (phylip/fasta) [format] \n \t --iter, -i : Phenotype shuffling iterations [iter] \n \t --phen, -p : Custom format phenotype file [phen] \n \t --out, -o : Outputfile redirective (not recommanded) [out]"
		sys.exit()

#################################################################################

# s -- SNP sequence file [seq]
# t -- Newick tree [tree]
# f -- SNP sequence file format (phylip/fasta) [format]
# i -- Phenotype shuffling iterations [iter]
# p -- Custom format phenotype file [phen]
# o -- Output file redirective (not recommanded) [out]

#################################################################################

ape=importr("ape")	# Required for phangorn
ph=importr("phangorn")	# Phylogenetic operations in R

print "All modules imported successfully"

t = PhyloTree(intree, alignment=seq, alg_format=seq_format)     # Main tree containing entire sequence
dtp = PhyloTree(intree)				 # Dummy tree for phenotype shuffling

print "Tree file read successfully"

phenfile=open(phen, "r")	# Phenotype file
phenlist=[]
for line in phenfile.readlines():
    phenlist.append([line.split("\t")[0].strip(), line.split("\t")[1].strip()])
phenfile.close()

phenotype={}	# Dictionary containing species names and their phenotype values

# Phenotype file can have three or four coloumns,
# Three column containing file : | SPECIES NAME | PHENOTYPE MEAN | STANDARD DEVIATION |
# Three column conataing file :  | SPECIES NAME | PHENOTYPE MEAN |


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



shphenotype=[]
shphenotype.append(phenotype.values())

#########################################################################################
# Reading complete ancestral sequence data generated through R (In the form of pseudo-patterns)
#########################################################################################
print "Reading complete ancestral sequence data generated through R"

rtree=PhyloTree(intree)	# Tree for "R" generated patterns
tree=ape.read_tree(intree)

rlist=[]
ropf = open("/tmp/rdata.dat", "r")	# rdata.dat is a rgp.R generated output file
for tab in ropf.readlines():
	tab=tab.rstrip()
	rlist.append(tab.split(" "))
ropf.close()

ori=np.array(rlist)

for node in rtree.traverse("postorder"):	# Patterns are being linked to their corresponding nodes
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
        node.data=map(lambda x: x, rlist[node.rtoken-1][1:])	# Its the sequence after name
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
print "Starting Phenotype manipulations"
                              
def phrange(node):
    ''' Recursive function to feed phenrange '''
    left, right = node.children
    if left.phenvalue is None:
        phrange(left)
    if right.phenvalue is None:
        phrange(right)

    if left.phenvalue < right.phenvalue :

        node.phenvalue=left.phenvalue + (dtp.get_distance(node, left)) * (right.phenvalue-left.phenvalue) / (dtp.get_distance(node, left) + dtp.get_distance(node, right))

    elif left.phenvalue > right.phenvalue :
        
        node.phenvalue=right.phenvalue + (dtp.get_distance(node, right)) * (left.phenvalue-right.phenvalue) / (dtp.get_distance(node, left) + dtp.get_distance(node, right))
    else :
        node.phenvalue=left.phenvalue

counter=0
for node in dtp.traverse("preorder"):
    node.add_features(counter=counter)
    counter += 1
    if node.is_leaf():
        node.add_features(phenvalue=phenotype[node.name][0])
    else:
        node.add_features(phenvalue=None) 

sflphen = []	# List to contain phenvalues for all the nodes after several shuffling attempts.

phrange(dtp)

for p in range(int(iterations)):
    shuffle(shphenotype[0])
    sflphen.append([])

    index = 0   
    for node in dtp.traverse("preorder"):     # serves two purpose, 1) Add existing values to sflphen 2) reset tree phenrange attribites
        sflphen[p].append(node.phenvalue)
        if node.is_leaf():
            node.phenvalue = shphenotype[0][index][0]
            index += 1 
        else:
            node.phenvalue = None
            
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
        else: pass
    else : pass


dt=rtree

pattern={}	# Contains all the possible patterns of SNPs across the tree, generated by rgp.R

sfldict={}	# Details of branches where nucleotide is changing and the corresponding p-value

bchanges={}	# Number of braches the nucleotide is changing and the corresponding pattern

#mafreq={}	# Stores minor allele frequency for individal SNP

counter=0
for node in dt.traverse('preorder'):
    node.add_features(counter=counter)
    counter += 1
        

for var in xrange(len(ref.data)): 
    currentdata=[]
    currentdict={"18" : 18}
    token=0
    for node in dt.traverse("postorder"):
        if node.is_leaf():
            try:
                currentdata.append(currentdict[node.data[var]])        # Pattern recognition
            except:
                currentdata.append(token)
                currentdict[node.data[var]]=token
                token+=1                
        else:
            pass   
        
    tdata=tuple(currentdata)

#    if tdata.count(0) < tdata.count(1):
#        af=tdata.count(0)*100/(tdata.count(0)+tdata.count(1))
#    elif tdata.count(0) > tdata.count(1):
#        af=tdata.count(1)*100/(tdata.count(0)+tdata.count(1))
#    else:
#        af=50

    try:
        pvalue=pattern[tdata]               
    except:
        temp=[]
        for node in dt.traverse('preorder'):       
            if not node.is_leaf():             
                left, right = node.children            
                if not left.data[var] == node.data[var]:         # Nucleotide substitution
                    if not left.data[var] == 18 :
                        temp.append([node.counter, left.counter])
                if not right.data[var] == node.data[var]:         # Nucleotide substitution
                    if not right.data[var] == 18 :
                        temp.append([node.counter, right.counter])
        
        temp.sort()
        tpp=map(tuple, temp)
        tp=tuple(tpp)
        dtemp=[]
        try:
            pvalue=sfldict[tp]
        except:
            for index in xrange(len(sflphen)):
                tsum=0.0
                if not len(temp) == 0:
                    for n in xrange(len(temp)):
                        tsum += abs(sflphen[index][temp[n][0]]-sflphen[index][temp[n][1]])
                    dtemp.append(math.pow((tsum/len(temp)),1))
                else :
                    dtemp.append('NA')
            if 'NA' in dtemp:
                pvalue=-1
            else:
                pvalue=scipy.stats.norm.sf((dtemp[0] - np.average(dtemp)) / np.std(dtemp))
            sfldict[tp]=pvalue
        pattern[tdata]=pvalue
        bchanges[tdata]=len(tp)
#        mafreq[tdata]=af
        
###############################################################################################
# Reading orignal sequence, evaluating by breaking into patterns
###############################################################################################
print "Reading orignal sequence, evaluating by breaking into patterns"


single=True

for node in t.traverse("postorder"):
    if single==True:
        if node.is_leaf():
            oriref=node
            single=False
        else: pass
    else : pass

#oriref=t.search_nodes(name="C57BL/6J")[0]

plist=[None for var in xrange(len(oriref.sequence))]
blist=[None for var in xrange(len(oriref.sequence))]
#maflist=[None for var in xrange(len(oriref.sequence))]

for var in xrange(len(oriref.sequence)):
    currentdata=[]
    currentdict={"-" : 18}
    token=0
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
              
    tdata=tuple(currentdata)
    
    try:
        plist[var]=pattern[tdata]
        blist[var]=bchanges[tdata]
#        maflist[var]=mafreq[tdata]
    except:
        print "Patterns missing"
###############################################################################################
try:
    opf=open(out, "w")
except:
    opf=open("/tmp/gpresult.dat", "w")

for i in range(len(plist)):
    opf.write(str(i)+'\t'+str(plist[i])+'\n')
opf.close()

print "Operation complete"
