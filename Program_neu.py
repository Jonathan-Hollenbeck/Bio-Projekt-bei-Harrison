import argparse
import os
import time
import datetime
import re
from difflib import SequenceMatcher

#declaring Methods

#read in all Arguments
def readInArguments():
    parser = argparse.ArgumentParser()

    #put in all the possible parameters
    parser.add_argument("--originalfasta", "-ogf",
        help = "The fasta file you wanna compare the others to.")
    parser.add_argument("--comparefastas", "-cpf", nargs = "+",
        help = "The fasta files you wanna compare to the original one.")
    parser.add_argument("--queries", "-q", nargs = "*", default = ["queries.txt"],
        help = "The genes you wanna compare. Either in the queries.txt or a list of specific genes without the >.")
    parser.add_argument("--threshold", "-t", default = "0",
        help="The threshold you wanna compare. Either a number with a percent symbole for percent calculation or a number for absolute calculation.")

    #make args global
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
    #looping through every line in the filestream
    for line in stream:
        #deleting line breaks
        line = line.replace("\n", "")
        #appending line
        fastaList.append(line)
    return fastaList

#parsing the content of a fasta in a dictionary
def parseToDict(stream):
    fileDict = {}
    keys = ""
    values = []
    items = ""
    #Reads the first input file line by line and distincts between >ID and Values
    #and writes the >IDs as keys and the Values as values in a dictionary
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

#checking if given string is a allowed fasta sequence
def isSequenceFasta(sequence):
    return bool(re.match("^[" + allowedFasta + "]+$", sequence))

#checking if file is in fasta format
def checkFastaFormat(file):
    #open filestream for given file  and paring into a list
    fileStream = openFileStream("input_files/" + file)
    fileContentAsList = parseFileToList(fileStream)
    fileStream.close()
    inGeneID = False
    inSequence = False
    #loop through all lines in the list
    for line in fileContentAsList:
        #ignore comment case
        if not line.startswith(";"):
            #checking for IDs
            if inGeneID == False and line.startswith(">"):
                inGeneID = True
                inSequence = False
            #checking if line after ID line is not an ID, but a sequence
            elif inGeneID == True and line.startswith(">"):
                print("ERROR: " + file + " is not in FASTA format, because there where 2 GenIDs in a row")
                exit(1)
            #if next line after ID is a sequence
            elif inGeneID == True:
                #checking if it is alphabetical
                if not isSequenceFasta(line):
                    print("ERROR: " + file + " is not in FASTA format, because of this line: " + line)
                    exit(1)
                inGeneID = False
                inSequence = True
            #if in sequence
            elif inSequence == True:
                #checking if it is alphabetical
                if not isSequenceFasta(line):
                    print("ERROR: " + file + " is not in FASTA format, because of this line: " + line)
                    exit(1)
            #if its none of the above print error
            else:
                print("ERROR: " + file + " is not in FASTA format, because of this line: " + line)
                exit(1)

# Checks the similarity of two given strings and returns it ratio
def similarRatio(seq1, seq2):
    return SequenceMatcher(None, seq1, seq2).ratio()

#Checks the difference of two given strings and returns it
def similarAbsolute(seq1, seq2):
    #loops through both sequences at once
    #and checks if the characters are the same
    #if not, add 1 to the difference
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

#comparing the fasta files using the queries
def compareFastasWithQueries(fasta1, fasta2, queries, threshold):
    #comparing via percentage
    if str(threshold).endswith("%"):
        #parsing threshold to float and without % sign
        threshold = float(threshold.replace("%", ""))
        #looping through every query and every ID in compare fasta
        for query in queries:
            for key in fasta2.keys():
                #if sequence length is equal, continue
                if len(fasta1.get(query)) == len(fasta2.get(key)):
                    #calculation the percentage similarity between the two sequences
                    similarity = similarRatio(fasta1.get(query), fasta2.get(key))
                    #checking if threshold was exceeded
                    if similarity*100 >= threshold:
                        #if was exceeded
                        #append the renamed ID with the sequence to all the outputs
                        appendToOutputs(query, key, str(similarity*100) + "%")
                        appendToFAOutput(query, fasta2.get(key))
                    else:
                        #if not
                        #append the original ID and sequence to all the output
                        appendToFAOutput(key, fasta2.get(key))
                else:
                    #if sequence was not equaly long
                    #append the original ID and sequence to the fasta output
                    appendToFAOutput(key, fasta2.get(key))
    #comparing via absolute value
    else:
        #parsing the threshold to a float
        threshold = float(threshold)
        #looping through every query and every ID in compare fasta
        for query in queries:
            for key in fasta2.keys():
                #if sequence is equaly long, continue
                if len(fasta1.get(query)) == len(fasta2.get(key)):
                    #calculate the absolue difference between the two sequences
                    difference = similarAbsolute(fasta1.get(query), fasta2.get(key))
                    #checking if threshold was exceeded
                    if difference <= threshold:
                        #if not
                        #append the renamed ID with the sequence to all the outputs
                        appendToOutputs(query, key, difference)
                        appendToFAOutput(query, fasta2.get(key))
                    else:
                        #if was exceeded
                        #append the original ID and sequence to all the output
                        appendToFAOutput(key, fasta2.get(key))
                else:
                    #if sequence was not equaly long
                    #append the original ID and sequence to the fasta output
                    appendToFAOutput(key, fasta2.get(key))

#manages the output lists
def appendToOutputs(query, key, simdiff):
    #append to fasta
    output_log.append("Changed ID " + key + " to " + query + "\n")
    #append to csv
    output_csv.append(key + "," + query + "," + str(simdiff) + "\n")

#manages the outputlist for the fasta output
def appendToFAOutput(query, sequence):
    #append query
    output_fa.append(query)
    #append sequence
    output_fa.append(sequence)

#creates a folder using the given directory,
#if it does not already exist
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('ERROR: Creating directory. ' + directory)

#writes all output files
def writeInOutputs(outputname):
    #check if output and logs folder are there and create it if not
    createFolder("./output/")
    createFolder("./logs/")
    #write the fasta file
    try:
        fo = open("./output/" + outputname.replace(".fa", "") + "_output.fa", "w+")
        #loop through all lines in the output list for the fasta
        for line in output_fa:
            #checks if line is a comment or a ID
            if line.startswith(">") or line.startswith(";"):
                #if yes just write it in
                fo.write(line + "\n")
            else:
                #if not break it in parts of length 60
                n = 60
                parts = [line[i:i + n] for i in range(0, len(line), n)]
                for part in parts:
                    fo.write(part + "\n")
        fo.close()
    except IOError:
        print("An error occured trying to write to " + outputname)

    #write in logfile
    #timestamp for logfilename
    currentTimestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H-%M-%S')
    try:
        lo = open("./logs/log_" + currentTimestamp + "_" + outputname.replace(".fa", "") + ".txt", "w+")
        lo.write("Program call was:\n")
        #write program call in log
        lo.write("python " + str(os.path.basename(__file__)) + str(args)[9:] + "\n")
        #loop through all lines in the output list for log and just write them in
        for element in output_log:
            lo.write(element)
        #write threshold and new outputfile in log
        lo.write("Chosen threshold was: " + str(args.threshold) + "\n" + "New filename for " + outputname + " is " + outputname.replace(".fa", "") + "_output.fa" + "\n")
        lo.close()
    except IOError:
        print("An error occured trying to append to log_" + currentTimestamp + ".txt")

    #write in csv file
    try:
        co = open("./output/" + outputname.replace(".fa", "") + "_csv.csv", "w+")
        #write similarity if percentage calculation was used
        if str(args.threshold).endswith("%"):
            co.write("oldID,newID,similarity\n")
        #write difference of absolute calculation was used
        else:
            co.write("oldID,newID,difference\n")
        #loop through all lines in the output list for csv and just write them in
        for element in output_csv:
            co.write(element)
        co.close()
    except IOError:
        print("An error occured trying to append to ./output/" + outputname.replace(".fa", "") + "_csv.csv")

#current time in milliseconds
current_milli_time = lambda: int(round(time.time() * 1000))

#actually doing stuff now

readInArguments();

#making a string with all allowed characters in a fasta sequence
#wich is every letter of the alphabet except J and O plus you can use * and -
allowedFasta = "AaBbCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtUuVvWwXxYyZz*-"

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
output_fa = []

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
    compareFastasWithQueries(ogFasta_dict, cpFasta_dict, args.queries, args.threshold)
    print("comparison between " + args.originalfasta + " and " + fasta + " completed!")

    print("writing output...")
    writeInOutputs(fasta)
    print("finished writing output!\n")

    output_csv = []
    output_log = []
    output_fa = []

    print("Time needed for this comparision: " + str((current_milli_time() - millisFasta)) + " milliseconds\n")

print("Time needed for all comparisions: " + str((current_milli_time() - millisAll)) + " milliseconds")
