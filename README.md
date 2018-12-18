# Finding and renaming existing ids of equal or similar sequences

Develop a tool that can be used to efficiently and accurately match transcript IDs between two or more transcriptomes.
The tool should take a list of transcripts, identify the same transcripts in another transcriptome and rename the transcripts so that they match the first set of IDs.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them
- Python version greater or equal than 3.6

Check your python version with:

```
python --version
```

### Usage

To use the program, simply change to project directory (where Program.py exists) with the console or terminal and
type

```
python Program_neu.py --help
```

to get the help message for the usage of this python program.

Some example program calls could be:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa
```

Here you have to put your query IDs in the queries.txt file inside the input_files folder.
Please put your Query IDs in there without the leading >.

You can query all IDs from first fasta file if you do a program call like:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa -q all
```

If you prefer checking just one ID i.e.: >ID1234 then you can do a program call like:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa -q ID1234
```

It is also possible to make a program call for more than one specific ID:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa -q ID1234 ID5678
```

If you would like to change the default threshold of 100% to i.e.: 98% then you can do a program call like:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa -t 98%
```

Or you can specify an absolute number of allowed differences:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa -t 2
```

means 2 nucleotides are allowed to be different.

You can do those checks for multiple fasta files too:

```
python Program_neu.py -ogf fasta_1.fa -cpf fasta_2.fa fasta_3.fa -t 2
```

All optional parameters can be used in the same program call.

### Resulting folder structure

Right after the program call, you will find the folders logs and output.

Inside the log folder is a resulting logfile with the timestamp of the program execution.
The logfile contains the whole program call with all parameters, every ID change, the chosen (or default threshold)
and the new fasta filename with the renamed sequences.

Inside the output folder you will find a file for each compared fasta file named with the old fasta filename and an appended _output.fa.
Additionally there is a corresponding _output_csv.csv with the old and new ID and their similarity.
