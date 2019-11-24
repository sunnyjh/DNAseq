from __future__ import division,print_function
import sys,os
import re
outfile = sys.argv[1]
fo = open(outfile,'w')
fo.write("chrom\tpos\tref\talt\ttotal\talt\tvaf\n")
bases = ["A","T","C","G","N","a","t","c","g","n"]
linecount = 0
for line in sys.stdin:
	linecount += 1
	if linecount % 3000 == 0:
		print("Processing "+str(linecount))
	tmp = line.strip().split("\t")
	# chr1    1138951 C       2503    .....	kkkkk	3328    .....	kkFFF
	total = 0
	mutations = {}
	deletions = {}
	insertions = {}
	for ba in "ATCGN":
		mutations[ba] = 0
	index = 3
	while index<len(tmp):
		total += int(tmp[index])	
		mapinfo=tmp[index +1]
		i = 0
		index += 3
		while i <len(mapinfo):
			mapbase=mapinfo[i]
			if mapbase=='.' or mapbase == "," or mapbase=="*":
				i += 1
			elif mapbase in bases:
				alt = mapbase.upper()
				mutations[alt] += 1
				i += 1
			elif mapbase=='$':
				i += 1
			elif mapbase=='^':
				if mapinfo[i+2] in bases:
					alt = mapinfo[i+2].upper()
					mutations[alt] += 1
				i += 3
			elif mapbase=='-':
				delnum=re.search('-(\d+)\w+',mapinfo[i:]).group(1)
				delalt=mapinfo[i+len(delnum)+1:i+len(delnum)+int(delnum)+1].upper()
				ref = tmp[2].upper() + delalt
				if ref not in deletions:
					deletions[ref] =1
				else:
					deletions[ref] += 1
				i=i+len(delnum)+int(delnum)+1
			elif mapbase=="+":
				insnum=re.search('\+(\d+)\w+',mapinfo[i:]).group(1)
				insalt=mapinfo[i+len(insnum)+1:i+len(insnum)+int(insnum)+1].upper()
				alt = tmp[2].upper() + insalt
				if alt not in insertions:
					insertions[alt] = 1
				else:
					insertions[alt] += 1
				i = i+len(insnum)+int(insnum)+1
			else:
				print(mapbase)
				i+=1
	flag = 1
	for base in mutations:
		alt_count = mutations[base]
		if alt_count==0:
			continue
		vaf = round(alt_count/total,6)
		olist =[str(x) for x in [tmp[0],tmp[1],tmp[2],base,total,alt_count,vaf]]
		fo.write("\t".join(olist)+"\n")
		flag = 0
	for ref in deletions:
		alt_count = deletions[ref]
		vaf = round(alt_count/total,6)
		# chr1    1138951 C
		olist = [str(x) for x in [tmp[0],tmp[1],ref,tmp[2],total,alt_count,vaf]]
		fo.write("\t".join(olist)+"\n")
		flag = 0
	for alt in insertions:
		alt_count = insertions[alt]
		vaf = round(alt_count/total,6)
		olist = [str(x) for x in [tmp[0], tmp[1], tmp[2],alt, total,alt_count,vaf]]
		fo.write("\t".join(olist)+"\n")
		flag = 0
	if flag:
		olist =[str(x) for x in [tmp[0],tmp[1],tmp[2],".",total,"0","0"]]
		fo.write("\t".join(olist)+"\n")
