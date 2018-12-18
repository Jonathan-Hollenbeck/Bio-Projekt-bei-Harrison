import argparse
import os
import time
import datetime
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
    parser.add_argument("--threshold", "-t", default = "0",
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

def appendToOutputs(query, key, simdiff):
    #append to fasta
    output_log.append("Changed ID " + key + " to " + query + "\n")
    #append to csv
    output_csv.append(key + "," + query + "," + str(simdiff) + "\n")

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('ERROR: Creating directory. ' + directory)

def writeInOutputs(fasta_dict, outputname):
    #check if output and logs folder are there and create it if not
    createFolder("./output/")
    createFolder("./logs/")
    #write the fasta file
    try:
        fo = open("./output/" + outputname.replace(".fa", "") + "_output.fa", "w+")
        for key in fasta_dict.keys():
            line = fasta_dict[key]
            n = 60
            newLine = [line[i:i + n] for i in range(0, len(line), n)]
            fo.write(key)
            fo.write("\n")
            for l in newLine:
                fo.write(l)
                fo.write("\n")
        fo.close()
    except IOError:
        print("An error occured trying to write to " + outputname)

    #write in logfile
    currentTimestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H-%M-%S')
    try:
        lo = open("./logs/log_" + currentTimestamp + "_" + outputname.replace(".fa", "") + ".txt", "w+")
        lo.write("Program call was:\n")
        lo.write("python " + str(os.path.basename(__file__)) + str(args)[9:] + "\n")
        for element in output_log:
            lo.write(element)
        lo.write("Chosen threshold was: " + str(args.threshold) + "\n" + "New filename for " + outputname + " is " + outputname.replace(".fa", "") + "_output.fa" + "\n")
        lo.close()
    except IOError:
        print("An error occured trying to append to log_" + currentTimestamp + ".txt")

    #write in csv file
    try:
        co = open("./output/" + outputname.replace(".fa", "") + "_csv.csv", "w+")
        if str(args.threshold).endswith("%"):
            co.write("oldID,newID,similarity\n")
        else:
            co.write("oldID,newID,difference\n")
        for element in output_csv:
            co.write(element)
        co.close()
    except IOError:
        print("An error occured trying to append to ./output/" + outputname.replace(".fa", "") + "_csv.csv")

def compareFastasWithQueries(fasta1, fasta2, queries, threshold):
    if str(threshold).endswith("%"):
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
                        appendToOutputs(query, key, str(similarity*100) + "%")
            else:
                difference = similarAbsolute(fasta1.get(query), fasta2.get(key))
                if difference <= threshold:
                    if not key in alreadyReplaced:
                        temp_dict[query] = temp_dict.pop(key)
                        alreadyReplaced.append(key)
                        appendToOutputs(query, key, difference)
    return temp_dict

#current time in milliseconds
current_milli_time = lambda: int(round(time.time() * 1000))

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
#error message if some queries are not in original file
for query in args.queries:
    if not query in ogFasta_dict.keys():
        print("ERROR: some or all queries are not in original File!")
        exit(1)
if not args.queries:
    print("ERROR: query list is empty! maybe check queries.txt.")
    exit(1)

#variable for log output and csv output
output_log = []
output_csv = []

millisAll = current_milli_time()

#looping through all compare fastas
for fasta in args.comparefastas:
    millisFasta = current_milli_time()
    checkFastaFormat(fasta)

    #reading in the compare Fasta file
    cpFastaStream = openFileStream("input_files/" + fasta)
    cpFasta_dict = parseToDict(cpFastaStream)
    cpFastaStream.close()

    print("starting to compare " + args.originalfasta + " and " + fasta + " with threshold " + str(args.threshold) + "...")
    writeInOutputs(compareFastasWithQueries(ogFasta_dict, cpFasta_dict, args.queries, args.threshold), fasta)
    print("comparison between " + args.originalfasta + " and " + fasta + " completed!")

    output_csv = []
    output_log = []

    print("Time needed for this comparision: " + str((current_milli_time() - millisFasta)) + " milliseconds\n")

print("Time needed for all comparisions: " + str((current_milli_time() - millisAll)) + " milliseconds")
