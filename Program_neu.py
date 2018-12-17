import argparse
from difflib import SequenceMatcher

#declaring Methods

#read in all Arguments
def readInArguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--originalfasta", "-ogf",
        help = "The fasta file you wanna compare the others to.")
    parser.add_argument("--comparefastas", "-cpf", nargs = "+",
        help = "The fasta files you wanna compare to the original one.")
    parser.add_argument("--queries", "-q", nargs = "*", default = ["queries.txt"],
        help = "The genes you wanna compare. Either in the queries.txt or a list of specific genes without the >.")
    parser.add_argument("--threshold", "-t", default = 0,
        help="The threshold you wanna compare. Either a number with a percent symbole for percent calculation or a number for absolute calculation.")

    global args

    args = parser.parse_args()

#opening a filestream for a given path
def openFileStream(file):
    try:
        return open(file, "r")
    except IOError:
        print("ERROR: An error occured trying to read " + file)
        exit(1)

#parsing the content of a file into a list
def parseFileToList(stream):
    fastaList = []
    for line in stream:
        line = line.replace("\n", "")
        fastaList.append(line)
    return fastaList

#parsing the content of a fasta in a dictionary
def parseToDict(stream):
    fileDict = {}
    keys = ""
    values = []
    items = ""
    # Reads the first input file line by line and distincts between >ID and Values
    # and writes the >IDs as keys and the Values as values in a dictionary
    for line in stream:
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
    return fileDict

#checking if file is in fasta format
def checkFastaFormat(file):
    fileStream = openFileStream("input_files/" + file)
    fileContentAsList = parseFileToList(fileStream)
    fileStream.close()
    inGeneID = False
    inSequence = False
    for line in fileContentAsList:
        if not line.startswith(";"):
            if inGeneID == False and line.startswith(">"):
                inGeneID = True
                inSequence = False
            elif inGeneID == True and line.startswith(">"):
                print("ERROR: " + file + "is not in FASTA format, because there where 2 GenIDs in a row")
                exit(1)
            elif inGeneID == True:
                if not line.isalpha():
                    print("ERROR: " + file + "is not in FASTA format, because of this line: " + line)
                    exit(1)
                inGeneID = False
                inSequence = True
            elif inSequence == True:
                if not line.isalpha():
                    print("ERROR: " + file + "is not in FASTA format, because of this line: " + line)
                    exit(1)
            else:
                print("ERROR: " + file + "is not in FASTA format, because of this line: " + line)
                exit(1)

# Checks the similarity of two given strings and returns it ratio
def similarRatio(seq1, seq2):
    return SequenceMatcher(None, seq1, seq2).ratio()

#Checks the difference of two given strings and returns it
def similarAbsolute(seq1, seq2):
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

#fasta[query] = fasta.pop(key)

def compareFastasWithQueries(fasta1, fasta2, queries, threshold):
    if threshold.endswith("%"):
        threshold = float(threshold.replace("%", ""))
        percent = True
    else:
        threshold = float(threshold)
        percent = False
    temp_dict = fasta2.copy()
    alreadyReplaced = []
    for query in queries:
        for key in fasta2.keys():
            if percent == True:
                similarity = similarRatio(fasta1.get(query), fasta2.get(key))
                if similarity*100 >= threshold:
                    if not key in alreadyReplaced:
                        temp_dict[query] = temp_dict.pop(key)
                        alreadyReplaced.append(key)
            else:
                difference = similarAbsolute(fasta1.get(query), fasta2.get(key))
                if difference <= threshold:
                    if not key in alreadyReplaced:
                        temp_dict[query] = temp_dict.pop(key)
                        alreadyReplaced.append(key)

#actually doing stuff now

readInArguments();

#checkingFastaFormat for original file
checkFastaFormat(args.originalfasta)

#reading in the original fasta file
ogFastaStream = openFileStream("input_files/" + args.originalfasta)
ogFasta_dict = parseToDict(ogFastaStream)
ogFastaStream.close()

#handling queries

#putting all genes from original fasta in queries
if "all" in args.queries:
    args.queries = []
    for key in ogFasta_dict.keys():
        args.queries.append(key.replace(">",""))
#if queries.txt was in -q, put all genes from it in queries but no doubles
elif "queries.txt" in args.queries:
    args.queries.remove("queries.txt")
    queriesStream = openFileStream("input_files/queries.txt")
    queries = parseFileToList(queriesStream)
    queriesStream.close()
    for query in queries:
        if query not in args.queries:
            args.queries.append(query)
#put > infront of every query
args.queries = [">" + query for query in args.queries]

#looping through all compare fastas
for fasta in args.comparefastas:
    checkFastaFormat(fasta)

    #reading in the compare Fasta file
    cpFastaStream = openFileStream("input_files/" + fasta)
    cpFasta_dict = parseToDict(cpFastaStream)
    cpFastaStream.close()

    print("starting to compare " + args.originalfasta + " and " + fasta + " with threshold " + args.threshold + "...")
    compareFastasWithQueries(ogFasta_dict, cpFasta_dict, args.queries, args.threshold)
    print("comparison between " + args.originalfasta + " and " + fasta + " completed!")
