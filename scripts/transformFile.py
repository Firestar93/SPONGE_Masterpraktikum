import numpy as np
import pandas as pd
import sys
import os

inDir = sys.argv[1]
mappingFile = sys.argv[2]
out = sys.argv[3]


def transformGenExpr(inPath, outPath):
    files = os.listdir(inPath)

    mapOfColumns = {}  # This will contain one list for each file, with the contents of column of expression values
    id = []
    fileCounter = 0

    for inFile in files:
        filename = inFile.split(".")[0]
        columnData = []
        f = open(inPath + inFile, 'r')
        first_line = f.readline()
        for line in f:
            b = line.split('\t')
            if fileCounter == 0:
                id.append(b[0])
            if len(b) >= 3:
                columnData.append(b[3].rstrip())
            else:
                # Row data is too short
                columnData.append('Error')
        f.close()
        mapOfColumns[filename] = columnData
        fileCounter = fileCounter + 1

    outFile = open(outPath, 'w')
    outFile.write("\t" + '\t'.join(id))
    outFile.write("\n")
    for key in mapOfColumns:
        outFile.write(key + "\t")
        outFile.write('\t'.join(mapOfColumns.get(key)))
        outFile.write('\n')
    outFile.close()

def transformMIRNAExpr(inPath, mappingPath, outPath):
    files = os.listdir(inPath)
    mapOfColumns = {}  # This will contain one list for each file, with the contents of column of expression values

    uniqId = getUniqMirna(inPath)
    mapping = getMapping(uniqId, mappingPath)

    for inFile in files:
        filename = inFile.split(".")[0].split("_")[1]
        columnData = {}
        f = open(inPath + inFile, 'r')
        for line in f:
            b = line.split('\t')
            if len(b) >= 2:
                columnData[b[0].rstrip()] = b[1].rstrip()
        f.close()

        result = []
        #get only needed data:
        for key in uniqId:
            if key in columnData.keys():
                result.append(columnData.get(key))
            else:
                result.append(0)
        mapOfColumns[filename] = result

    outFile = open(outPath, 'w')
    outFile.write("\t" + '\t'.join(mapping.get(key) for key in uniqId))
    outFile.write("\n")
    for key in mapOfColumns:
        outFile.write(key + "\t")
        outFile.write('\t'.join(str(x) for x in mapOfColumns.get(key)))
        outFile.write('\n')
    outFile.close()

def getUniqMirna(inPath):
    files = os.listdir(inPath)
    uniqId = {}
    for inFile in files:
        f = open(inPath + inFile, "r")
        for line in f:
            b = line.split("\t")
            if b[0].rstrip() in uniqId.keys():
                uniqId[b[0].rstrip()] = uniqId.get(b[0].rstrip()) + 1
            else:
                uniqId[b[0].rstrip()] = 1
        f.close()

    #filter mirna, that is <20 times in the data found
    return [k for k, v in uniqId.items() if v > 20]

def getMapping(uniqMirna, mappingPath):
    mappingMap = {}
    f = open(mappingPath, "r")
    for line in f:
        b = line.split("\t")
        mimat = b[0].rstrip()

        for id in b[1:-1]:
            mappingMap[id.rstrip()] = mimat

    return mappingMap

    f.close()
########################################################################################################################
#Methodenaufrufe
########################################################################################################################
transformGenExpr(inDir, out)
#transformMIRNAExpr(inDir, mappingFile,out)





