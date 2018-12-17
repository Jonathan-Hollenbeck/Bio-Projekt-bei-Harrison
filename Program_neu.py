import argparse

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

def checkFastaFormat(fileContentAsList, filename):
    inGeneID = False
    inSequence = False
    for line in fileContentAsList:
        if not line.startswith(";"):
            if inGeneID == False and line.startswith(">"):
                inGeneID = True
                inSequence = False
            elif inGeneID == True and line.startswith(">"):
                print("ERROR: " + filename + "is not in FASTA format, because there where 2 GenIDs in a row")
                exit(1)
            elif inGeneID == True:
                if not line.isalpha():
                    print("ERROR: " + filename + "is not in FASTA format, because of this line: " + line)
                    exit(1)
                inGeneID = False
                inSequence = True
            elif inSequence == True:
                if not line.isalpha():
                    print("ERROR: " + filename + "is not in FASTA format, because of this line: " + line)
                    exit(1)
            else:
                print("ERROR: " + filename + "is not in FASTA format, because of this line: " + line)
                exit(1)


def compareFastas(fasta1ContentAsList, fasta2ContentAsList):
    print("compareFastas")

#actually doing stuff now

readInArguments();

#reading in the original fasta file
ogFastaStream = openFileStream("input_files/" + args.originalfasta)
ogFasta = parseFileToList(ogFastaStream)
ogFastaStream.close()

checkFastaFormat(ogFasta, "input_files/" + args.originalfasta)

#putting all genes from originalfasta in queries
if "all" in args.queries:
    args.queries = []
    for line in ogFasta:
        if line.startswith(">"):
            args.queries.append(line)
#if queries.txt was in -q, put all genes from it in queries but no doubles
elif "queries.txt" in args.queries:
    args.queries.remove("queries.txt")
    queriesStream = openFileStream("input_files/queries.txt")
    queries = parseFileToList(queriesStream)
    queriesStream.close()
    for query in queries:
        if query not in args.queries:
            args.queries.append(query)

#looping throu all compare fastas
for fasta in args.comparefastas:
    cpFastaStream = openFileStream("input_files/" + fasta)
    cpFasta = parseFileToList(cpFastaStream)
    cpFastaStream.close()

    checkFastaFormat(cpFasta, "input_files/" + fasta)
