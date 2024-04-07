import sys

outputH = open(sys.argv[2], "w")

to_write = dict()
for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    to_write[line[0] + line[1][1:] + line[2]] = "\t".join((line[0], line[1][1:], line[2] + "\n"))
    to_write[line[0] + line[3] + line[4]] = "\t".join((line[0], line[3], line[4] + "\n"))

for key in to_write:
    outputH.write(to_write[key])
