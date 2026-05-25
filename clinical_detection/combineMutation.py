from __future__ import print_function, division
import sys, os


def parseVCF(infile):
    records = {}
    sample = os.path.basename(infile).split(".")[0]
    with open(infile) as fi:
        for line in fi:
            tmp = line.strip().split("\t")
            if line.startswith("#"):
                continue
            else:
                #chr1    11856378        Mut1    G       A       .       .       SP=0.000908;ST=1.28     GT:AD:SI:DP     0/1:516,407:287,229,181,226:923
                k = ":".join([tmp[0], tmp[1], tmp[3], tmp[4]])
                records[k] = 1
    return records


def sortkey(key):
    info = key.split(":")
    chr2int = {}
    for i in range(23):
        chr2int["chr" + str(i)] = i
    chr2int["chrX"] = 23
    chr2int["chrY"] = 24
    chr2int['chrM'] = 25
    if info[0] not in chr2int:
        chr2int[info[0]] = 26
    return [chr2int[info[0]], int(info[1])]


def parseMutation(fi, sample, keeps):
    tick = 0
    mutations = {}
    depth = {}
    total_muts = {}
    #if sample not in mutations:
    #	mutations = {}
    #	depth = {}
    # chr1    1138951 C       .       903     0       0
    previos_pos = current_pos = ""
    while 1:
        file_num = fi.tell()
        line = fi.readline()
        tmp = line.strip().split("\t")
        #file_num = fi.tell()
        if tmp[0] == "chrom": continue
        if len(tmp) < 3:
            print("error:"+sample,tmp)
            break
        tmp[2] = tmp[2].upper()
        tmp[3] = tmp[3].upper()
        key = ":".join(tmp[0:2])
        #print(keeps)
        if key not in keeps:
            #print(sample,key,file=sys.stderr)
            fi.seek(file_num)
            break
        current_pos = key
        if key not in depth:
            depth[key] = int(tmp[4])
        if tmp[2] == ".":
            pass
        else:
            key = ":".join(tmp[0:4])
            mutations[key] = [tmp[5], tmp[4], tmp[6]]
            total_muts[key] = 1
        previos_pos = key
    #else:
    #    print('s_pos',current_pos,file=sys.stderr)
    return mutations, depth, total_muts

    for line in fi:
        #file_num = fi.tell()
        tmp = line.strip().split("\t")
        if tmp[0] == "chrom":
            continue
        key = ":".join(tmp[0:2])

        current_pos = key
        if current_pos != previos_pos:
            tick += 1
        previos_pos = key
        #print(sample,previos_pos,current_pos)
        if key not in depth:
            depth[key] = int(tmp[4])
        if tmp[2] == ".":
            pass
        else:
            key = ":".join(tmp[0:4])
            mutations[key] = [tmp[5], tmp[4], tmp[6]]
            total_muts[key] = 1
        if tick >= offset:
            #fi.seek(file_num)
            break
    return mutations, depth, total_muts, tick >= offset


def load_bed(bedfile):
    regions = []
    for line in open(bedfile):
        tmp = line.strip().split("\t")
        #chrM    1464    1584
        regions.append(tmp)
    return regions

def main():
    infile, bedfile, outfile = sys.argv[1:]
    regions = load_bed(bedfile)
    chr2int = {}
    for i in range(23):
        chr2int["chr" + str(i)] = i
    chr2int["chrX"] = 23
    chr2int["chrY"] = 24
    chr2int['chrM'] = 25
    mutations = {}
    depth = {}
    total_muts = {}

    sample2record = {}
    sample2fi = {}
    samples = []
    for line in open(infile):
        tmp = line.strip().split(" ")
        # YF1901P YF1901P.bgm.cons.xls
        sample2fi[tmp[0]] = open(tmp[1], 'r')
        #parseMutation(tmp[1],tmp[0])
        sample2record[tmp[0]] = parseVCF(tmp[2])
        samples.append(tmp[0])
    fo = open(outfile, 'w')
    fo.write("chrom\tpos\tref\talt\ttotal_count\talt_count\tvaf\t" +
             "\t".join(samples) + "\n")
    flag = 1
    keeps = {}
    for r in regions:
        for i in range(int(r[1])+1,int(r[2])+1):
            s = ":".join([r[0],str(i)])
            keeps[s] = 1
        mutations = {}  #update
        depth = {}
        total_muts = {}
        flags = []
        
        for sample in sample2fi:
            mutations[sample], depth[sample], t_muts = parseMutation(
                sample2fi[sample], sample, keeps = keeps)
            #flags.append(f)
            for k in t_muts:
                total_muts[k] = 1
            #print(sample,depth[sample])
            #print(sample,"mutations", mutations[sample])
            #print(sample,"depth", depth[sample])
        for k in sorted(total_muts, key=lambda x: sortkey(x)):
            total_count = alt_count = vaf = 0
            olist = []
            newkey = ":".join(k.split(":")[0:2])
            #print(k, newkey)
            for sample in samples:
                if k in mutations[sample]:
                    infos = mutations[sample][k]
                    total_count += int(infos[1])
                    alt_current = 0
                    if k in sample2record[sample]:
                        alt_count += 0
                        alt_current = 0
                        vaf = 0
                    else:
                        alt_current = int(infos[0])
                        alt_count += int(infos[0])
                        try:
                            vaf = round(alt_current / int(infos[1]), 7)
                        except:
                            vaf = 0
                    olist.append("%s/%s,%s" % (infos[0], infos[1], vaf))
                elif newkey in depth[sample]:
                    total_count += depth[sample][newkey]
                    olist.append("0/%s,0" % (depth[sample][newkey]))
                else:
                    #print(k, sample)
                    raise (KeyError(k + "," + sample + " not found"))

            if total_count >= 1:
                vaf = round(alt_count / total_count, 7)
            else:
                vaf = 0
            olist = ["\t".join(k.split(":")), total_count, alt_count, vaf
                     ] + olist
            olist = [str(x) for x in olist]
            fo.write("\t".join(olist) + "\n")

if __name__ == "__main__":
    main()
