import sys

input_h = open(sys.argv[1])
output_h = open(sys.argv[2], "w")
# input_h.readline() # skip header
for line in input_h:
    if not line.startswith("#"):
        line = line.strip().split()
        try:
            start = line[9]
        except IndexError:
            print(line)
            raise IndexError

        stop = line[10]
        strand = line[8]
        if strand == "-":
            tmp = stop
            stop = start
            start = tmp
        output_h.write("\t".join((line[2], start, stop, line[0], ".", strand)) + "\n") # inverted line[0] and [2] bc its te seq THEN chr

# chr1    10001   10468   (TAACCC)n       28.10659867662235       +
