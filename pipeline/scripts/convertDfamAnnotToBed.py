import sys

input_h = open(sys.argv[1])
output_h = open(sys.argv[2], "w")
input_h.readline() # skip header
for line in input_h:
    line = line.strip().split("\t")
    start = line[9]
    stop = line[10]
    strand = line[8]
    if strand == "-":
        tmp = stop
        stop = start
        start = tmp
    output_h.write("\t".join((line[0], start, stop, line[2], ".", strand)) + "\n")

# chr1    10001   10468   (TAACCC)n       28.10659867662235       +