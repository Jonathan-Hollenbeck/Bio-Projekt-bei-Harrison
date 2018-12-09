import time
import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fasta1")
parser.add_argument("fasta2")
parser.add_argument("queries")
parser.add_argument("output")
args = parser.parse_args()

fasta1 = "input_files/" + args.fasta1
fasta2 = "input_files/" + args.fasta2
queries = "input_files/" + args.queries
output = "output/" + args.output

try:
    f1 = open(fasta1, "r")
except IOError:
    print("An error occured trying to read file 1")
try:
    f2 = open(fasta2, "r")
except IOError:
    print("An error occured trying to read file 2")
try:
    queries = open(queries, "r")
except IOError:
    print("An error occured trying to read queries.txt")

# Declaring and initialising some dictionaries and key variables
f1_dict = {}
f2_dict = {}
keys1 = ""
keys2 = ""

def parseToDict(file):
    fileDict = {}
    keys = ""
    # Reads the first input file line by line and distincts between >ID and Values
    # and writes the >IDs as keys and the Values as values in a dictionary
    for line in file:
        if line.startswith(">"):
            values = []
            items = line.split("\n")
            keys = items[0]
            fileDict[keys] = values
        else:
            items = line.split("\n")
            values.append(items[0])
            str = ''.join(values)
            fileDict[keys] = str
    file.close()
    return fileDict

def checkAndRenameID(dict1, dict2, id):
    # rewrite ID of same sequences in the second input file
    f1_dict = dict1
    f2_dict = dict2
    inputID = id
    temp_dict = {}

    # check if dictionary 1 has more or equal items than dictionary 2
    if len(f1_dict.values()) >= len(f2_dict.values()):
        # iterate through all items of dictionary 2 (file2)
        for y in list(f2_dict):
            if f1_dict[inputID] == f2_dict[y]:
                print(inputID)
                print(y)
                temp_dict.update({inputID: f2_dict[y]})
                logToOutput(y, inputID)
            else:
                temp_dict.update({y: f2_dict[y]})
    else:
        for x in list(f2_dict):
            if f2_dict[x] == f1_dict[inputID]:
                print(inputID)
                print(x)
                temp_dict.update({inputID: f2_dict[x]})
                logToOutput(x, inputID)
            else:
                temp_dict.update({x: f2_dict[x]})
    return temp_dict

def writeToOutput(dict):
    try:
        lo = open("output/log.txt", "a")
    except IOError:
        print("An error occured trying to append to log.txt")
    lo.write("New filename for " + fasta2 + " is " +  output + "\n")
    lo.close()
    try:
        fo = open(output, "w")
    except IOError:
        print("An error occured trying to write to file")
    for x in dict:
        line = dict[x]
        n = 60
        newLine = [line[i:i + n] for i in range(0, len(line), n)]
        fo.write(x)
        fo.write("\n")
        for l in newLine:
            fo.write(l)
            fo.write("\n")
    fo.close()

def logToOutput(oldId, newId):
    try:
        lo = open("output/log.txt", "a")
    except IOError:
        print("An error occured trying to append to log.txt")
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    lo.write(st + ":\n" + "\tChanged ID\t" + oldId + "\n\t\tto\t" + newId + "\n")
    lo.close()

f1_dict = parseToDict(f1)
f2_dict = parseToDict(f2)

queryIDs = parseToDict(queries).keys()

for ID in queryIDs:
    f2_dict = checkAndRenameID(f1_dict, f2_dict, ID)

writeToOutput(f2_dict)