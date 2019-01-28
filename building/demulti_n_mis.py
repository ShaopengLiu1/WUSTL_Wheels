# for now count statistics only
# regard mismach, no deletion / insertion allowed
import pandas as pd
import time
import sys

read1=sys.argv[1]
read2=sys.argv[2]
index=sys.argv[3]
mis_cutoff=int(sys.argv[4])

#barcode=sys.argv[5]  map the number into dict and replace str to dict key

# find number of match for SAME LENGTH string
def dist_xy (seq1, seq2):
    mis=len(seq1)-sum([ a==b for (a,b) in zip(seq1, seq2)])
    return mis


def find_match (seq, ref, mismatch):
    mis_record=[]
    for start in range(len(ref)-len(seq)+1):
        mis_record.append(dist_xy(seq, ref[start:start+len(seq)]))
    if min(mis_record) <= mismatch:
        return 1
    else:
        return 0

w=open(index, 'r')
H3=0
FLAG=0
FatC_FLAG=0
GP=0
FatC_Trl=0
Trl=0
unmatch=0

i=0
while i <= 1e20:
    i1 = w.readline()
    i2 = w.readline().strip('\n')
    i3 = w.readline()
    i4 = w.readline()
    if i1 == '': break
    if find_match("CAGATC", i2, mis_cutoff):
        Trl += 1
    elif find_match("CTTGTA", i2, mis_cutoff):
        GP += 1
    elif find_match("CGATGT", i2, mis_cutoff):
        FatC_Trl += 1
    elif find_match("ACAGTG", i2, mis_cutoff):
        H3 += 1
    elif find_match("GCCAAT", i2, mis_cutoff):
        FLAG += 1
    elif find_match("GTGAAACG", i2, mis_cutoff):
        FatC_FLAG += 1 
    else: 
        unmatch += 1
    i += 4
w.close()

d={'H3': [H3], 'FLAG':[FLAG], 'FatC_FLAG':[FatC_FLAG], 'GP':[GP], 'FatC_Trl':[FatC_Trl], 'Trl':[Trl], 'unmatch':[unmatch]  }
out=pd.DataFrame(d).transpose()
out.columns=['count']
out['ratio']=out['count']*100 / out['count'].sum()
out['item']=list(out.index.values)
out=out[['item','count','ratio']]
out.to_csv("demultiplex_record_%s_mismatch.txt" %mis_cutoff, sep="\t", index=False)





