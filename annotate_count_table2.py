#!/usr/bin/env python
import sys
from collections import defaultdict

def make_ensid_list(filename):
    ensid_dict=defaultdict(list)
    annot_dict=defaultdict(list)
    for line in open(filename):
        l= line.split("\t")
        #print l
        gene_id=l[0].rstrip()
        #print gene_id 
        gene_name=l[1].rstrip()
        if len(gene_name)<1:
            gene_name='NA'
        #print gene_name
        gene_description=l[2].rstrip()
        if len(gene_description) <1 :
            gene_description='NA'
        #print gene_description
        #print l
        ensid_dict[gene_id]=gene_name
        annot_dict[gene_id].append(gene_description)
        #print gene_id, ensid_dict[gene_id], annot_dict[gene_id]
    return ensid_dict, annot_dict

def go_thru_count_table(filename, ensid_dict, annot_dict):
    ct_nid = 0 
    out=open('annotated_'+filename, "w")
    for line in open(filename):
        if (not line.startswith("NM") and not line.startswith("ENS")):
            out.write(line)
        #elif line.startswith('NM')or line.startswith('ENS'):
        if line.startswith('ENS'):
            #print line, len(line)
            l=line.rstrip().split('\t')
            #print len(l)
            #print l[0]
            if l[0].rstrip() in ensid_dict: 
                #print 'yes in dict'
                #print l[0]
                if l[0].rstrip() in annot_dict:
                    an=set(annot_dict[l[0]])
                    an2 = list(an)
                    #print an2
                    out.write(ensid_dict[l[0]].rstrip()+'\t'+'\t'.join(an2)+'\t'+line)
                else:
                    out.write(ensid_dict[l[0]].rstrip()+'\t'+line)
            #out.write('\n') d
            else: 
                #print 'not in dict', l[0]
                ct_nid+=1
    print 'not in dict=', ct_nid 
        
if __name__ == "__main__":
    (ensid_dict, annot_dict)=make_ensid_list(sys.argv[1])
    go_thru_count_table(sys.argv[2], ensid_dict, annot_dict)
