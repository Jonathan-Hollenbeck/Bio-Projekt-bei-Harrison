import time
import datetime
import argparse
from difflib import SequenceMatcher

#constructor for making the sequenz of code clear
def constructor():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1")
    parser.add_argument("fasta2")
    parser.add_argument("queries")
    parser.add_argument("output")
    parser.add_argument("threshold")
    args = parser.parse_args()

    #referenz all global variables used in this function
    global fasta1
    global fasta2
    global queries
    global queryIDs
    global output
    global threshold
    global f1copy
    global f2copy
    global f1
    global f2
    global f1_dict
    global f2_dict
    global keys1
    global keys2
    global ts
    global st

    fasta1 = "input_files/" + args.fasta1
    fasta2 = "input_files/" + args.fasta2

    #manage queries arg
    if args.queries[-4:] == ".txt":
        queries = "input_files/" + args.queries
        queries = openFile(queries)
        queryIDs = parseToDict(queries).keys()
    else:
        queryIDs = {args.queries: ""}

    #manage output arg
    if args.output == "":
        output = "output/" + args.fasta2[:-3] + "_output.fa"
    else:
        output = args.output
    threshold = float(args.threshold)


    #open needed filestreams
    f1copy = openFile(fasta1)
    f2copy = openFile(fasta2)
    f1 = openFile(fasta1)
    f2 = openFile(fasta2)

    # Declaring and initialising some dictionaries and key variables
    f1_dict = {}
    f2_dict = {}
    keys1 = ""
    keys2 = ""

    # Timestamp for execution
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H-%M-%S')

#function to open a filestream. returns filestream
def openFile(file):
    try:
        return open(file, "r")
    except IOError:
        print("An error occured trying to read " + file)

#checks if the given file is a fasta or not
def checkFastaFormat():
    inGeneID = 0
    inSequence = 0
    for line in f1copy:
        justText = line.split("\n")
        if inGeneID == 0 and justText[0].startswith(">"):
            inGeneID = 1
        elif inGeneID == 1:
            if not justText[0].isalpha():
                print("Error occured: Given Sequence contains other symbols")
                exit(1)
            inGeneID = 0
            inSequence = 1
        elif inSequence == 1 and justText[0].startswith(">"):
            inSequence = 0
            inGeneID = 1
        elif inSequence == 1:
            inSequence = 1
            if not justText[0].isalpha():
                print("Error occured: Given Sequence contains other symbols")
                exit(1)
        else:
            print("Error occured: Given file is not in fasta format.")
            exit(1)

    try:
        f2copy = open(fasta2, "r")
    except IOError:
        print("An error occured trying to read file 1")
    for line in f2copy:
        justText = line.split("\n")
        if inGeneID == 0 and justText[0].startswith(">"):
            inGeneID = 1
        elif inGeneID == 1:
            if not justText[0].isalpha():
                print("Error occured: Given Sequence contains other symbols")
                exit(1)
            inGeneID = 0
            inSequence = 1
        elif inSequence == 1 and justText[0].startswith(">"):
            inSequence = 0
            inGeneID = 1
        elif inSequence == 1:
            inSequence = 1
            if not justText[0].isalpha():
                print("Error occured: Given Sequence contains other symbols")
                exit(1)
        else:
            print("Error occured: Given file is not in fasta format.")
            exit(1)

# Checks the similarity of two given strings and returns it ratio
def similar(seq1, seq2):
    return SequenceMatcher(None, seq1, seq2).ratio()

#parses the content from a file into a dictionary
def parseToDict(file):
    fileDict = {}
    keys = ""
    values = []
    items = ""
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

# rewrite ID of same sequences in the second input file
def checkAndRenameID(dict1, dict2, id):
    f1_dict = dict1
    f2_dict = dict2
    inputID = id
    temp_dict = {}

    # check if dictionary 1 has more or equal items than dictionary 2
    if len(f1_dict.values()) >= len(f2_dict.values()):
        # iterate through all items of dictionary 2 (file2)
        for y in list(f2_dict):
            # saves the similarity between sequence1 and sequence2
            similarity = similar(f1_dict[inputID], f2_dict[y])
            if f1_dict[inputID] == f2_dict[y]:
                print(inputID)
                print(y)
                temp_dict.update({inputID: f2_dict[y]})
                logToOutput(y, inputID)
            elif similarity >= threshold:
                print(inputID)
                print(y)
                print(threshold)
                temp_dict.update({inputID: f2_dict[y]})
                logToOutput(y, inputID + " has similarity of: " + str(similarity))
            else:
                temp_dict.update({y: f2_dict[y]})
    else:
        for x in list(f2_dict):
            # saves the similarity between sequence1 and sequence2
            similarity = similar(f1_dict[inputID], f2_dict[x])
            if f2_dict[x] == f1_dict[inputID]:
                print(inputID)
                print(x)
                temp_dict.update({inputID: f2_dict[x]})
                logToOutput(x, inputID)
            elif similarity >= threshold:
                print(inputID)
                print(x)
                print(threshold)
                temp_dict.update({inputID: f2_dict[x]})
                logToOutput(x, inputID + " has similarity of: " + str(similarity))
            else:
                temp_dict.update({x: f2_dict[x]})
    return temp_dict

def writeToOutput(dict):
    #write in logfile
    try:
        lo = open("logs/log_" + st + ".txt", "w+")
        lo.write("Chosen threshold was: " + str(threshold) + "\n" + "New filename for " + fasta2 + " is " + output + "\n")
        lo.close()
    except IOError:
        print("An error occured trying to append to log_" + st + ".txt")

    #write in outputfile
    try:
        fo = open(output, "w")
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
    except IOError:
        print("An error occured trying to write to " + output)

#writes changes that were made into the logfile
def logToOutput(oldId, newId):
    try:
        lo = open("logs/log_" + st + ".txt", "a")
    except IOError:
        print("An error occured trying to append to log_" + st + ".txt")
    lo.write("\n" + "\tChanged ID\t" + oldId + "\n\t\tto\t" + newId + "\n")
    lo.close()

constructor()

checkFastaFormat()

f1_dict = parseToDict(f1)
f2_dict = parseToDict(f2)

for ID in queryIDs:
    f2_dict = checkAndRenameID(f1_dict, f2_dict, ID)

writeToOutput(f2_dict)
