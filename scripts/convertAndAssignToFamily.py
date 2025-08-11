# arg 1 is ucsc input
# arg 2 is bed output
# arg 3 is family table output

import sys

output = open(sys.argv[2], "w")
familiesTE = open(sys.argv[3], "w")
file_content = list()
maxScore = 0
families = dict() # associate element ID to big family

#bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand repName repClass repFamily repStart repEnd repLeft id 
inputHandler = open(sys.argv[1])

inputHandler.readline()  # skip first line

for line in inputHandler:
    line = line.strip().split()
    if maxScore < int(float(line[1])):  # score is the second column
        maxScore = int(float(line[1]))
    if line[10] not in families:  # repName is the key
        families[line[10]] = [line[11] + '/' + line[12], 0, 0]
        # initialize counters + merge repClass_repFamily
    families[line[10]][1] += 1
    families[line[10]][2] += int(line[7]) - int(line[6])
    if line[9] == "C":  # strand is in 9th position
        line[9] = "-"

    # seqname / start / end /
    line = [line[5], line[6], line[7], line[10], line[1], line[9]]
    file_content.append(line)

for index in range(0, len(file_content)):
    file_content[index][4] = \
        str(float(file_content[index][4]) / maxScore * 1000)
    file_content[index] = "\t".join(file_content[index]) + "\n"

output.writelines(file_content)

for value in families:
    familiesTE.write(value +
                     "\t" +
                     families[value][0] +
                     "\t" +
                     str(families[value][1]) +
                     "\t" +
                     str(families[value][2] / families[value][1]) + "\n")
