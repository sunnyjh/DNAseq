from __future__ import print_function,division
import sys,os
infile, outfile = sys.argv[1:]
#YF1901P YF1901P.bgm.cons.xls
fo = open(outfile,'w')
chr2int = {}
for i in range(23):
	chr2int["chr"+str(i)] = i
chr2int["chrX"] = 23
chr2int["chrY"] = 24

mutations = {}
depth = {}
total_muts = {}

def sortkey(key):
	info = key.split(":")
	return [chr2int[info[0]],int(info[1])]
	
def parseMutation(infi,sample):
	if sample not in mutations:
		mutations[sample] = {}
		depth[sample] = {}
	# chr1    1138951 C       .       903     0       0
	for line in open(infi):
		tmp = line.strip().split("\t")
		if tmp[0] == "chrom":
			continue
		key = ":".join(tmp[0:2])
		if key not in depth[sample]:
			depth[sample][key] = int(tmp[4])
		if tmp[2] == ".":
			pass
		else:
			key = ":".join(tmp[0:4])
			mutations[sample][key] = [tmp[5],tmp[4],tmp[6]]
			total_muts[key] = 1

for line in open(infile):
	tmp = line.strip().split(" ")
	# YF1901P YF1901P.bgm.cons.xls
	parseMutation(tmp[1],tmp[0])

samples = mutations.keys()
fo.write("chrom\tpos\tref\talt\ttotal_count\talt_count\tvaf\t"+"\t".join(samples)+"\n")
for k in sorted(total_muts,key= lambda x:sortkey(x)):
	total_count = alt_count = vaf = 0
	olist = []
	newkey = ":".join(k.split(":")[0:2])
	for sample in samples:
		if k in mutations[sample]:
			infos = mutations[sample][k]
			total_count += int(infos[1])
			alt_count += int(infos[0])
			if int(infos[1])>=1:
				vaf = round(int(infos[0])/int(infos[1]),7)
			else:
				vaf = 0
			olist.append("%s/%s,%s" % (infos[0],infos[1],vaf))
		elif newkey in depth[sample]:
			total_count += depth[sample][newkey]
			olist.append("0/%s,0" % (depth[sample][newkey]))
		else:
			print(k,sample)
			raise
	if total_count>=1:
		vaf = round(alt_count/total_count,7)
	else:
		vaf = 0
	olist = ["\t".join(k.split(":")),total_count,alt_count,vaf] + olist
	olist = [str(x) for x in olist]
	fo.write("\t".join(olist)+"\n")

