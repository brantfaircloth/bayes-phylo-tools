import pdb 

counts = {}

for line in open('917_loci/917_sub_models.txt', 'rU'):
    v = line.strip().split('\t')[1]
    print v
    #pdb.set_trace()
    counts.setdefault(v,0)
    counts[v] += 1
for k,v in counts.iteritems():
    print k,v
