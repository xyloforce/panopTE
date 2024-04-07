import sys
#bin    swScore milliDiv        milliDel        milliIns        genoName        genoStart       genoEnd genoLeft        strand  repName repClass        repFamily       repStart        repEnd  repLeft id
input_h = open(sys.argv[1])
input_h.readline()
input_h.readline()
input_h.readline() # skip the first three lines
output = open(sys.argv[2], "w")
output.write("\t".join(("#bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id")) + "\n")
count = 0
for line in input_h:
    line = line.split()
    line[7] = line[7][1:-1]
    line.insert(0, str(count))
    output.write("\t".join(line) + "\n")
    count += 1