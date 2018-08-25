#!/usr/bin/env python

import sys
import itertools
from collections import defaultdict


def combinations_by_subset(seq, r):
	if r:
		for i in xrange(r - 1, len(seq)):
			for cl in combinations_by_subset(seq[:i], r - 1):
				yield cl + (seq[i],)
	else:
		yield tuple()

def rename_list(names):
	#rename to R names of headers
	new_names=[]
	for n in names: 
		N=n.replace('+','.').replace('-','.')
		#print N
		new_names.append(N)
	return new_names


def do_all(ctrl,expt,fileName,RNames,orderLineDict1, orderLineDict2, number_of_samples):
	# Case 1 : 2x2 expt vs control 
	expt_comb=combinations_by_subset(RNames, int(number_of_samples))
	for combo in expt_comb: 
		#print "combo:",combo
		# load count table
		print "ct <- read.table('"+fileName+"', h=T, row.names=1)"
		#remove unecessary columns
		cols_removed=[]
		for Name in RNames: 
			#print orderLineDict2[combo]
			if Name not in combo: 
				print "ct$"+Name+" <- NULL"
				cols_removed.append(Name)
		#print 'removed samples=',cols_removed
		# make design
		NamesAfterRem= [] 
		for name in RNames:
			if not name in cols_removed :
				NamesAfterRem.append(name)
		#print 'Names after removal', NamesAfterRem
		factor=""
		for i in range (0,len(NamesAfterRem)):
			if NamesAfterRem[i] in RexptNames:
				if i == (len(NamesAfterRem)-1):
					factor +='"expt"'
				else : 
					factor +='"expt",'
			if NamesAfterRem[i] in RctrlNames: 
				if i == (len(NamesAfterRem)-1):
					factor +='"ctrl"'
				else : 
					factor +='"ctrl",'
		#print factor
		design_names='_'.join(map(str, list(combo)))
		#print design_names

		print 'design.'+design_names+'=data.frame(row.names=colnames(ct), type=factor(c('+factor+')))'

		# make deseq2 object
		print 'diffexp.dsq =  DESeqDataSetFromMatrix(countData=ct, colData=design.'+design_names+', design=~type)'
		print 'diffexp.dsq = DESeq(diffexp.dsq)'
		print 'diffexp.res = results(diffexp.dsq, contrast=c("type", "expt","ctrl"))'
		print 'diffexp.res_up = diffexp.res[diffexp.res$log2FoldChange > 0 & !is.na(diffexp.res$padj) & diffexp.res$padj < 0.05,]'
		print 'diffexp.res_dn = diffexp.res[diffexp.res$log2FoldChange < 0 & !is.na(diffexp.res$padj) & diffexp.res$padj < 0.05,]'
		print 'write.table(as.data.frame(diffexp.res), file="DifferentialExpression_'+design_names+'.txt",sep="\\t", quote=FALSE)'
		print 'write.table(as.data.frame(diffexp.res_up), file="DifferentialExpression_'+design_names+'up.txt",sep="\\t", quote=FALSE)'
		print 'write.table(as.data.frame(diffexp.res_dn), file="DifferentialExpression_'+design_names+'dn.txt",sep="\\t", quote=FALSE)'

if __name__ == "__main__":

	# arguments to give 
	Names=sys.argv[1].split(',')
	#print 'All the headers:', Names
	exptNames= sys.argv[2].split(',')		# experimental conditions
	#print 'Names of experiments:', exptNames, len(exptNames)
	ctrlNames= sys.argv[3].split(',')			# control conditions
	#print 'Names of controls:', ctrlNames, len(ctrlNames), '\n\n'
	fileName = str(sys.argv[4])					# name of file

	RNames=rename_list(Names)
	RexptNames=rename_list(exptNames)
	RctrlNames=rename_list(ctrlNames)

	#print Names
	ct=0
	expt=[]
	ctrl=[]
	orderLineDict1=defaultdict(str)
	orderLineDict2=defaultdict(str)
	for RName in RNames: 
		ct+=1
		orderLineDict1[ct]=RName
		orderLineDict2[RName]=ct
		#print ct, Name
		if RName in RexptNames: 
			expt.append(ct)
			#print RName, ct, 'expt!'
		elif RName in RctrlNames: 
			ctrl.append(ct)
			#print RName, ct, 'ctrl!'
	#print 'controls:', ctrl
	#print 'expts', expt

	
	for i in (3, len(Names)-1):
		#print i,"!!!"
		do_all(ctrl,expt,fileName,RNames,orderLineDict1, orderLineDict2, i)
	
