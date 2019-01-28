# for now count statistics only
import pandas as pd
import sys
import random
import time

read1=sys.argv[1]
read2=sys.argv[2]
index=sys.argv[3]   # the index file
barcode=sys.argv[4] # the barcode match, should contain ONLY 2 column: sample_name barcode

def save_dict(dict, name="demultiplex_record_0_mis.txt"):
    out=pd.DataFrame(dict, index=[0]).transpose()
    out.columns=['count']
    out['ratio']=out['count']*100 / out['count'].sum()
    out['item']=list(out.index.values)
    out=out[['item','count','ratio']]
    out.to_csv(name, sep="\t", index=False)


match_table=pd.read_table(barcode, header=None, names=['sample','barcode'])
barcode_list=list(match_table.loc[:, 'barcode'])
item_list=list(match_table.loc[:,'sample'])
dict_match={key:value for (key, value) in zip(item_list, barcode_list)}

tS = time.time()
# count distri and write to each output files
with open(index, 'r') as w:
    with open(read1,'r') as r1:
        with open(read2, 'r') as r2:
            dict_count={key:0 for key in item_list} # intialize start count
            i=0
            unmatch=0
            f={}
            for num, sample in enumerate(item_list):
                f[str(num+1)]=file('%s_1.fastq' %sample,'w')  # the number is 0-indexed
                f[str(num*-1)]=file('%s_2.fastq' %sample,'w')   # to distinguish str2
            f['unmatch_1']=file('unmatch_1.fastq','w')
            f['unmatch_2']=file('unmatch_2.fastq','w')
            f['unfound_index']=file('unfound_index.fastq','w')
            while i <= 1e20:
                try:
                    i1 = w.readline()
                    i2 = w.readline()
                    i3 = w.readline()
                    i4 = w.readline()
                    L1 = r1.readline()
                    L2 = r1.readline()
                    L3 = r1.readline()
                    L4 = r1.readline()
                    R1 = r2.readline()
                    R2 = r2.readline()
                    R3 = r2.readline()
                    R4 = r2.readline()
                    if i1 == '': break
                    i+=4
                    bingo=0
                    for num, sample in enumerate(item_list):
                        if dict_match[sample] in i2:
                            dict_count[sample]+=1
                            bingo=1
                            f[str(num+1)].write(L1+L2+L3+L4)
                            f[str(num*-1)].write(R1+R2+R3+R4)
                            break
                    if bingo==0:
                        unmatch+=1
                        f['unmatch_1'].write(L1+L2+L3+L4)
                        f['unmatch_2'].write(R1+R2+R3+R4)
                        f['unfound_index'].write(i1+i2+i3+i4)
                    dict_count['unmatch']=unmatch
                except:
                    print 'All done'
                    break
            for key in f:
                f[str(key)].close()
save_dict(dict_count, "demultiplex_record_0_mis.txt")

# random order input seq to confirm there is no ambiguity
random.shuffle(item_list)
with open(index, 'r') as w:
    dict_count={key:0 for key in item_list} # intialize start count
    i=0
    unmatch=0
    while i <= 1e20:
        try:
            i1 = w.readline()
            i2 = w.readline()
            i3 = w.readline()
            i4 = w.readline()
            if i1 == '': break
            i+=4
            unmatch+=1
            for sample in item_list:
                if dict_match[sample] in i2:
                    dict_count[sample]+=1
                    unmatch+=-1
                    break
            dict_count['unmatch']=unmatch
        except:
            print 'All done'
            break
save_dict(dict_count, "random_order_demultiplex_record_0_mis.txt")


tE = time.time()
print 'Cost ',(tE-tS),' sec'

