from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from itertools import chain, compress
import pandas as pd
import numpy as np
import re


def zeroForUnknown(x):
    match = re.search('Unknown*', x)
    if(match):
        x = "Un: Unknown"
    return x


# getting the index from scafMOA_matrix
def getMatrixIndex(scaffoldList, compoundSMILES):
    matchInd = []
    cmpdMol = Chem.MolFromSmiles(compoundSMILES)
    for index in range(len(scaffoldList)):
        scaffold = Chem.MolFromSmiles(scaffoldList[index])
        if(cmpdMol.HasSubstructMatch(scaffold)):
            matchInd.append(index)     
    if(len(matchInd) > 0):

        return matchInd
    else:

        return 'NA'


def subStructureSearch(aSMILE, SMILESlist):
    
    indexList = []
    for index in range(len(SMILESlist)):
        subStr = Chem.MolFromSmiles(SMILESlist[index])
        if subStr.HasSubstructMatch(Chem.MolFromSmiles(aSMILE)):
            indexList.append(index)

    return indexList


def getScaffold_MOA_matrix(cmpdsSMILES, cmpdsMOA, scaffoldSMILES):
    moa = []
    matchIndex = []
    scaffold_moa_col = []
    cmpdList = cmpdsSMILES
    scaffolds = scaffoldSMILES
    cmpdsMOAlist = cmpdsMOA
   
    for aScaffold in scaffolds:
        indexList = subStructureSearch(aScaffold, cmpdList)
        matchIndex.append(indexList)
        moa.append(pd.unique(cmpdsMOAlist[indexList]))
    
    moaList = pd.Series(pd.unique(cmpdsMOAlist))
    moaList = moaList.apply(zeroForUnknown)
    moaList = pd.unique((sorted(moaList)))
    scfId = ["SCF" + str(number) for number in list(range(1, len(scaffolds) + 1))]

    for moaType in moaList:
        scaffold_moa_col.append(moaType.split(":")[0])
    
    zeroMatrix = np.zeros((len(scfId) * len(moaList),), dtype=np.int).reshape(len(scfId), len(moaList))
    dfa = pd.DataFrame({'SCF_ID': scfId, 'SMILES': scaffolds})#order problem later
    dfb = pd.DataFrame(zeroMatrix, columns = scaffold_moa_col)
    scafMOA_matrix = pd.concat([dfa, dfb], axis = 1)
    #all possibilities based on scaffold_substructure_MOA mapping
    for row in range(len(moa)):
        for amode in moa[row]:
            amode = zeroForUnknown(amode)
            col = amode.split(":")[0]
            scafMOA_matrix.loc[row, col] = 1

    return scafMOA_matrix



def getMatrixIndex(scaffolds_From_scafMOA_matrix, compoundSMILES):
    matchInd = []
    scaffoldList = scaffolds_From_scafMOA_matrix
    cmpdMol = Chem.MolFromSmiles(compoundSMILES)
    for index in range(len(scaffoldList)):
        scaffold = Chem.MolFromSmiles(scaffoldList[index])
        if(cmpdMol.HasSubstructMatch(scaffold)):
            matchInd.append(index)     
    if(len(matchInd) > 0):

        return matchInd
    else:

        return 'NA'



def getMOA(cmpdSMILES, scafMOA_matrix):

    moaList = []
    moa = ['NA']
    scaffolds = scafMOA_matrix['SMILES']
    matrixInd = getMatrixIndex(scaffolds, cmpdSMILES)
    if(matrixInd != 'NA'):
        selectedMatrix = scafMOA_matrix.iloc[matrixInd , 2:]
        for ind in matrixInd:
            moatypes = selectedMatrix.loc[ind, :][selectedMatrix.loc[ind, :] == 1] 
            moaList.append(moatypes.index.values.tolist())
        moaList = list(chain.from_iterable(moaList))
        freqList = [moaList.count(aMOA) for aMOA in set(moaList)]
        moa = list(compress(set(moaList),[frequency == max(freqList) for frequency in freqList]))
        UnremovedMOA = list(compress(moa, [aMOA != 'Un' for aMOA in moa]))
        if(UnremovedMOA != []):
            moa = UnremovedMOA
        else: moa = ['NA']
    else: moa = ['NA']
    
    return moa


def getMOAlist(molecule, scafMOA_matrix):

    moaList = []
    df = pd.DataFrame()
    for acmpd in molecule:
        moa = getMOA(acmpd, scafMOA_matrix)
        if(len(moa) > 1):
            newstr = moa[0]
            for ind in range(1, len(moa)):
                newstr += ',' + moa[ind]
            moa = newstr
        else:
            moa = moa[0]        
        moaList.append(moa)
    df['SMILES']= molecule
    df['MOA']= moaList

    return df