import sys

output = open(sys.argv[2], "w")
for line in open(sys.argv[1]):
    line = line.strip().split()
    new_id = "_".join((line[0], line[1], line[2]))
    size = int(line[2]) - int(line[1])
    # we need :
    # 50 bp after end
    # at last 100 inside
    if not size >= 150:
        if line[5] == "+":
            line[1] = int(line[1]) - (150 - size) # ensure there is enough space on the border
            line[2] = int(line[2]) + 50 # go outside
        else:
            line[2] = int(line[2]) + (150 - size)
            line[1] = int(line[1]) - 50
    else: # already the right size
        if line[5] == "+":
            line[2] = int(line[2]) + 50 # go outside
        else:
            line[1] = int(line[1]) - 50
    
    line[3] = new_id
    if int(line[1]) < 0:
        line[1] = 0
    line[1] = str(line[1])
    line[2] = str(line[2])
    output.write("\t".join(line) + "\n")