# p4_bygene_compotest.py v1.0
#
# Chris Laumer 3 April 2014 GPL v2
#
# This is a minimally written script, designed to individually test a set of genes 
# for compositional heterogeneity using null simulations as implemented in p4 v0.9.0, 
# simulating on inferred ML gene trees. Run it in a folder for which you have, with 
# each gene, the following information:
#
# a.) A phylip-formatted MSA, named [gene].phylip
# b.) A maximum likelihood tree in newick form, here called RAxML_bestTree.[gene]_ML
# c.) A model, in the format output by ProteinModelSelection.pl, named [gene]_model.txt.
#     The only models supported are those used by p4; you may have to modify ProteinModelSelection.pl
#     accordingly.
#
# If these conditions are met, running the script should print the name of each gene, plus
# the p-value for the null test, both to standard output and to a file called "by_gene_compo_test.txt".
# If you are running it on a computer with many CPUs, the python multiprocessing module
# will automatically try to use all of those processes (good for an HPC cluster, bad for a desktop).


import glob
from p4 import *
from multiprocessing import Pool, current_process

genes = glob.glob('OG*.phylip')

for i in range(len(genes)):
	genes[i] = genes[i].rstrip('.phylip')

def main(gene):
	var.doCheckForDuplicateSequences = False
	read(gene + '.phylip')
	read('RAxML_bestTree.' + gene + '_ML')
	a = var.alignments[0]
	t = var.trees[0]
	
	d = Data()
	
	t.data = d
	
	handle = open(gene + '_model.txt')
	model = handle.read().rstrip('\n')
	handle.close()
	
	def map_to_p4_modelnames(str):
		'''
		Takes in a string of the RAxML model name output by ProteinModelSelection.pl, outputs a tuple with the corresponding p4 name, and a T/F for whether to use empirical frequencies
		'''
		mod_map = [('DAYHOFF','d78'), ('LG', 'lg'), ('WAG', 'wag'), ('JTT', 'jtt'), ('MTREV', 'mtREV24'), ('RTREV', 'rtRev'), ('MTMAM', 'mtmam'), ('CPREV','cpREV'), ('RTREV','rtRev'), ('BLOSUM62','blosum62'), ('HIVB','hivb'), ('MTART','mtart'), ('MTZOA','mtzoa')] 
		if str[-1] == 'F' and str != 'DAYHOFF':
			basemod = str[:len(str)-1]
			emp = True
		else:
			basemod = str
			emp = False
		for mod in mod_map:
			if basemod == mod[0]:
				p4_mod = mod[1]
		return (p4_mod,emp)
	
	# Set model parameters
	p4_model = map_to_p4_modelnames(model)
	if p4_model[1]:
		t.newComp(free=1, spec='empirical')
		t.newRMatrix(free=0, spec=p4_model[0])
	else:
		t.newComp(free=0, spec=p4_model[0])
		t.newRMatrix(free=0, spec=p4_model[0])
	t.setNGammaCat(nGammaCat=4)
	t.newGdasrv(free=1, val=0.5)
	t.setPInvar(free=0, val=0.0)
	
	t.optLogLike(verbose=0)
	
	#Avoid stochastic failures...
	status = False
	while status is False: 
		try:
			ret = t.compoTestUsingSimulations(doIndividualSequences=0,nSims=100,verbose=0, doChiSquare=0)	
			print gene, ret
			status = True
		except:
			status = False
	var.alignments = []
	var.trees = []
	handle = open('by_gene_compo_test.txt', 'a+')
	handle.write(gene + '\t' + str(ret) + '\n')
	handle.close()

if __name__ == '__main__':
	pool = Pool()
	calc = pool.imap(main, genes)
	pool.close()
	pool.join()	
