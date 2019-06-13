#!/usr/bin/python
from __future__ import division
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import reverse_complement
import pandas
from collections import defaultdict
import subprocess
import sys
import shutil


def main():
    options = process_arguments()

    # Import the sample data from the information in the config file
    samples = process_config_file(options.input)

    # Pairwise comparison of each sample against each other
    # 1. split the input fasta files
    # 2. Make blast databases of the split input fasta files
    # 3. Pairwise blast the split intput fasta files
    sliceDB = sliceAndCreateDB(samples, options)

    # Take each slice in each sample and blast against all the other samples. If "no hit" is detected
    # against other sample databases, then keep it as "unique". If hits are detected, then determine
    # if the slice hits only members of particular groupings and assign appropriately.
    if (options.noblast != True):
        blastresults = blastSlices(samples, options) #this should be empty, gets filled in next section
    else:
        blastresults = defaultdict(lambda: defaultdict(dict))


    blastresults, pairwiseresults = processPairwiseBlastResults(samples, options, blastresults, sliceDB)


    # Filtering Method checkpoint - sorts out the datatype issue with the dataframe.
    results = pandas.read_csv("pairwisedata."+str(options.size)+"."+str(options.identity)+"."+str(options.minlength)+".csv", index_col=['slice'],dtype=str)

    final = finalOutputProcessing(results, options, samples, blastresults, sliceDB)
    print final


'''
SampleData Class
'''
class SampleData(object):
    '''
    Describes the sample data.
    This could be one or more sequences that constitute a sample and/or genome
    Genbank or Fasta could be redundant as required
    '''
    def __init__(self, name, fastafilename=None, groups=None):
        self.name = name
        self.fastafilename = fastafilename
        self.groups = groups
        self.sequences = self.readSequences(fastafilename)

    def readSequences(self, fastafilename):
        return [seq for seq in SeqIO.parse(fastafilename, "fasta")]

    def __str__(self):
        outline = ['Sample data: %s' % self.name]
        outline.append('Fastafile: %s' % self.fastafilename)
        outline.append('Groups: %s' % self.groups)
        outline.append('Num Sequences: %s' % len(self.sequences))
        return os.linesep.join(outline) + os.linesep

'''
GENERAL FUNCTIONS
'''
def finalOutputProcessing(results, options, samples, blastresults, sliceDB):
    try:
        os.mkdir("./" + options.output)
    except:
        shutil.rmtree("./" + options.output, ignore_errors=True)
        os.mkdir("./" + options.output)
    for columnName, columnData in results.iteritems():
        if 'GroupVariable' in columnName:
            unique = columnData.unique()
            variableFolder = "./" + options.output + "/" + str(columnName)
            try:
                os.mkdir(variableFolder)
            except:
                pass
            for variable in unique:
                print("Processing %s %s" %(columnName,variable))
                try:
                    os.mkdir(variableFolder + "/" + str(variable))
                except:
                    pass
                positiveSampleList, negativeSampleList = getSampleLists(columnName, variable, samples)
                if (len(negativeSampleList) > 0):
                    positiveDF = results.ix[:,positiveSampleList]
                    negativeDF = results.ix[:,negativeSampleList]
                    for posData in positiveDF.itertuples():
                        sliceName = posData[0]
                        if 'Negative' not in posData:
                            negData = negativeDF.loc[sliceName]
                            negDataList = negData.tolist()
                            if 'positive' not in negDataList:
                                sequence = sliceDB[sliceName]["sequence"]
                                sequence.description = sliceName
                                associatedSequences = getAssociatedSequences(sequence, sliceName, blastresults, samples)
                                numHits = len(associatedSequences) - 1
                                sliceFastaOutfile = open(variableFolder + "/" + str(variable) + "/"+ sliceName + "." + str(numHits) + "hits.fasta", "w")
                                SeqIO.write(associatedSequences, sliceFastaOutfile, "fasta")
                else:
                    positiveDF = results.ix[:,positiveSampleList]
                    for posData in positiveDF.itertuples():
                        sliceName = posData[0]
                        if 'Negative' not in posData:
                            sequence = sliceDB[sliceName]["sequence"]
                            sequence.description = sliceName
                            associatedSequences = getAssociatedSequences(sequence, sliceName, blastresults, samples)
                            numHits = len(associatedSequences) - 1
                            sliceFastaOutfile = open(variableFolder + "/" + str(variable) + "/" + sliceName + "." + str(numHits) + "hits.fasta", "w")
                            SeqIO.write(associatedSequences, sliceFastaOutfile, "fasta")
    return "FINISHED!"


def getAssociatedSequences(seq, sliceName, blastresults, samples):
    associatedSequences = []
    associatedSequences.append(seq)
    seqblastresults = blastresults[sliceName]
    for hit in seqblastresults:
        hitstart = int(seqblastresults[hit]["hitstart"])
        hitend = int(seqblastresults[hit]["hitend"])
        for sample in samples:
            for sequence in sample.sequences:
                if sequence.id == hit:
                    if (hitstart > hitend):
                        start = hitend
                        end = hitstart
                        associatedSlice = sequence[start:end]
                        thisid = sequence.id
                        associatedSlice = associatedSlice.reverse_complement()
                        associatedSlice.id = sample.name + "_" + thisid
                        associatedSlice.description = "[" + str(start) + "-" + str(end) + "]"
                    else:
                        start = hitstart
                        end = hitend
                        associatedSlice = sequence[start:end]
                        associatedSlice.description = "[" + str(start) + "-" + str(end) + "]"
                        thisid = sequence.id
                        associatedSlice.id = sample.name + "_" + thisid
                    associatedSequences.append(associatedSlice)
    return associatedSequences


def processPairwiseBlastResults(samples, options, blastresults, sliceDB):
    size = int(options.size)
    pairwiseresults = pandas.DataFrame()
    for sample in samples:
        print("Processing sample %s" % sample.name)
        for line in [l.strip() for l in open("blasts/slices." + str(size) + "." + sample.name + ".megablast", "rU")]:
            query, hit, pid, alignmentlength, mismatches, gaps, qstart, qend, hitstart, hitend, evalue, bitscore = line.split(
                '\t')
            if (int(alignmentlength) >= int(options.minlength)):
                if (float(pid) >= float(options.identity)):
                    blastresults[query][hit]["hit"] = hit
                    blastresults[query][hit]["pid"] = pid
                    blastresults[query][hit]["alignmentlength"] = alignmentlength
                    blastresults[query][hit]["mismatches"] = mismatches
                    blastresults[query][hit]["gaps"] = gaps
                    blastresults[query][hit]["qstart"] = qstart
                    blastresults[query][hit]["qend"] = qend
                    blastresults[query][hit]["hitstart"] = hitstart
                    blastresults[query][hit]["hitend"] = hitend
                    blastresults[query][hit]["evalue"] = evalue
                    blastresults[query][hit]["bitscore"] = bitscore
                    blastresults[query][hit]["groups"] = sliceDB[query]["groups"]
                    blastresults[query][hit]["qsample"] = sliceDB[query]["sample"]
                    blastresults[query][hit]["qsequence"] = sliceDB[query]["sequence"]
                    i = 1
                    for group in sliceDB[query]["groups"]:
                        pairwiseresults.set_value(query, "GroupVariable" + str(i), group)
                        i += 1
                    pairwiseresults.set_value(query, sample.name, "positive")
                    pairwiseresults.set_value(query, "slice", query)
    pairwiseresults = pairwiseresults.set_index('slice')
    pairwiseresults.to_csv(
        "pairwisedata." + str(size) + "." + str(options.identity) + "." + str(options.minlength) + ".csv",
        na_rep="Negative")
    return [blastresults, pairwiseresults]


def blastSlices(samples, options):
    size = int(options.size)
    blastresults = defaultdict(lambda: defaultdict(dict))
    try:
        os.mkdir("./blasts")
    except:
        shutil.rmtree("./blasts", ignore_errors=True)
        os.mkdir("./blasts")
    for sample in samples:
        print("Searching sample %s" % sample.name)
        blastcmd = "blastn -query slices." + str(size) + ".fasta -db databases/" + sample.name + " -out blasts/slices." + str(size) + "." + sample.name + ".megablast -outfmt 6 -evalue 0.001 -max_target_seqs 1000 -num_threads " + str(options.threads)
        pairwisechild = subprocess.Popen(str(blastcmd),
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             universal_newlines=True,
                                             shell=(sys.platform != "win32"))
        pairwiseoutput, pairwiseerror = pairwisechild.communicate()
    return blastresults


def sliceAndCreateDB(samples, options):
    size = int(options.size)
    sliceDB = defaultdict(lambda: defaultdict(dict))
    outseqs = open("slices." + str(size) + ".fasta", "w")
    count = 0
    for sample in samples:
        for seq in sample.sequences:
            for i in range(0, len(seq), size):
                slice = seq[i:i + size]
                slice.id = "slice" + str(count)
                SeqIO.write(slice, outseqs, "fasta")
                sliceDB[slice.id]["sample"] = sample.name
                sliceDB[slice.id]["sequence"] = seq[i:i + size]
                sliceDB[slice.id]["sequence_start"] = i
                sliceDB[slice.id]["sequence_stop"] = i + size
                sliceDB[slice.id]["groups"] = sample.groups
                count += 1
        if (options.noblast != True):
            blastdbcmd = "makeblastdb -in " + sample.fastafilename + " -out databases/" + sample.name + " -dbtype nucl"
            blastdbchild = subprocess.Popen(str(blastdbcmd),
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE,
                                            universal_newlines=True,
                                            shell=(sys.platform != "win32"))
            blastdboutput, blastdberror = blastdbchild.communicate()
    return sliceDB


def getSampleLists(groupColumnName, groupVariable, samples):
    columnlist = groupColumnName.split('GroupVariable')
    column = int(columnlist[1]) - 1
    positiveList = []
    negativeList = []
    for sample in samples:
        groups = sample.groups
        if groups[column] == groupVariable:
            positiveList.append(sample.name)
        else:
            negativeList.append(sample.name)
    return positiveList, negativeList


def process_config_file(configfile):
    samples = []
    for line in [l.strip() for l in open(configfile, 'rU') if l.strip() and not l.startswith('#')]:
            sample, groups, fastafile = re.split('\t', line)
            grouplist = groups.split(',')
            samples.append(SampleData(sample, fastafile, grouplist))
    return samples


def process_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Config file describing sequence files and groupings', required=True)
    parser.add_argument('--size', help='Desired size (in nucleotides)', required=True, type=int, default=1000)
    parser.add_argument('--identity', help='Minimum percentage identity for splitting positive/negative sequences [0-100]', required=False, default=95.0, type=float)
    parser.add_argument('--output', help='Output folder name for the results', required=False)
    parser.add_argument('--threads', help='Number of threads available', type=int, required=False, default=2)
    parser.add_argument('--minlength', help='Minimum blast hit length (PERCENTAGE) to consider a MATCH to another genome. Good hits with length less than this will be described as unique between the samples. Default=0.5', type=int, default=0.5, required=True)
    parser.add_argument('--noblast', help='For debugging', action="store_true")
    return parser.parse_args()

if __name__ == '__main__':
    main()