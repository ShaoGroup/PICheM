#!/usr/bin/python
# Extract Compounds' Features by RDKit.
# Songpeng Zu /zusongpeng@gmail.com/

#-- Import packages

# Import numpy and pandas for data structures in python
import numpy as np
from pandas import DataFrame, Series
import pandas as pd

# Import rdkit package for chemical feature extraction.
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
# Import 2D Pharmacophore Fingerprints
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
# Import Fingerprints.
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs, Torsions
# Import Descriptor Calculation
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
# Import urllib2 to get the the compounds from CHEMBL directly.
import urllib2
## Import the regular expression in python.
#import re

#-- Define a global variable to restrict the lines to read at one time.
MAXLINES = 1000

#-- Definitions of functions.
def readFILE(f):
    with open(filen) as f:
        compounds = f.readlines()
    f.close()
    return(compounds)

# This function doesnot work well, it seems that the chemilid and the comid are not the same as recorded in the website.
def getSMILESfromCHEMBL(chemblid):
    """Generate the SMILES format of a compound recorded as CHEMBL ID.
    The record is like "CHEMBLE116"

    Keyword arguments:
    chemblid: the compound's CHEMBL ID.

    Return: a string represented the SMILES format.
    """
    webstr = "https://www.ebi.ac.uk/chembldb/download_helper/getmol/"
    comid = chemblid.replace("CHEMBL","")
    cominfor = urllib2.urlopen(webstr+comid)
    d = cominfor.read()
    m = Chem.MolFromMolBlock(cominfor.read())
    return(Chem.MolToSmiles(m))

def phychemDataFrame(chempandas,namecol,smicol):
    """ Generate the physicochemical properties of the compounds.
    The compounds are stored in the DataFrame Structure defined by pandas.

    Keyword arguments:
    chempandas: the compounds stored in DataFrame, which contain the name and SMILES as columns.
    namecol: the column number of the name of SMILES.
    smicol: the column number of SMILES in the DataFrame.

    Return: a DataFrame of the compounds merging chempadas and the phychem by columns. 
    If None is detected given a SMILES-like string, it would be not deleted. 
    Note: The SMILES output by Chem.MolToSmiles is canonical, and might be different with the original.
    Add the names to different compounds.
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
#   molsmi = [x for x in molsmitmp if x is not None]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1

    nms = [x[0] for x in Descriptors._descList]
    nms.remove("MolWt")
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    descrs = [calc.CalcDescriptors(x) for x in molsmi]
    descrsmat = np.matrix(descrs)
    df = DataFrame(descrsmat, index = [x.GetProp("_Name") for x in molsmi], columns=nms)
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)

def phychemFileStream(inputfile,namecol,smicol,outputfile,chunknum=1000):
    """ Generate the physicochemical properties of the compounds from files.
    The compounds are recorded as the SMILES format.

    Keyword argument:
    inputfile: the file of compounds in which each row a compound.
    namecol: the column number of the name of SMILES.
    smicol: the column number of SMILES in the DataFrame.
    outputfile: the file of results in which each row a compound and each column a physiclchemical feature.
    chunknum: number of lines are read once a time due to the file might be quite large.
              default is 1000, should be smaller than MAXLINES.
    
    Return: outputfile containing these properties.
    """
    assert chunknum <= MAXLINES
    chunker = pd.read_table(inputfile,sep="\t",chunksize=chunknum)
    nms = [x[0] for x in Descriptors._descList]
    nms.remove("MolWt")
#    title = "\t".join(nms.append("SMILES"))+"\n"
#    f = open(outputfile,"w")
#    f.write(title)
#    f.close()
    for piece in chunker:
        df = phychemDataFrame(piece,namecol,smicol)
#       dfnew = df.append(Series(list(df.columns.values),index=list(df.columns.values)),ignore_index=True)
        df.to_csv(outputfile,mode="a",header=False,sep="\t",index=False)
    return(0)

def MACCfpDataFrame(chempandas,namecol,smicol):
    """ Generate the physicochemical properties of the compounds.
    The compounds are stored in the DataFrame Structure defined by pandas.

    Keyword arguments:
    chempandas: the compounds stored in DataFrame, which contain the name and SMILES as columns.
    namecol: the column number of the name of SMILES.
    smicol: the column number of SMILES in the DataFrame.

    Return: a DataFrame of the compounds merging chempadas and the fingerprints by columns. 
    If None is detected given a SMILES-like string, it would be not deleted. 
    Note: The SMILES output by Chem.MolToSmiles is canonical, and might be different with the original.
    Add the names to different compounds.
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1
    # MACC Fingerprints.
    fps = [MACCSkeys.GenMACCSKeys(x) for x in molsmi]
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) # how to name the col?
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)

def TOPOLOGYfpDataFrame(chempandas,namecol,smicol):
    """
    Topology-based fingerprints 2048 bits. 
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1
    # TOPOLOGY Fingerprints.
    fps = [FingerprintMols.FingerprintMol(x) for x in molsmi]
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) # how to name the col?
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)
    
def ATOMPAIRSfpDataFrame(chempandas,namecol,smicol):
    """
    AtomPairs-based fingerprints 2048 bits.
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1
    # ATOMPAIRS Fingerprints.
    fps = [Pairs.GetAtomPairFingerprintAsBitVect(x) for x in molsmi]
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) # how to name the col?
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)

def TORSIONSfpDataFrame(chempandas,namecol,smicol):
    """
    Torsions-based fingerprints 2048 bits. 
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1
    # TORSIONS Fingerprints.
    fps = [Torsions.GetTopologicalTorsionFingerprintAsIntVect(x) for x in molsmi]
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) # how to name the col?
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)

def MORGANfpDataFrame(chempandas,namecol,smicol):
    """
    MORGAN-based fingerprints 1024 bits. 
    """
    assert chempandas.shape[0] <= MAXLINES
    molsmitmp = [Chem.MolFromSmiles(x) for x in chempandas.iloc[:,smicol]]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chempandas.iloc[i,namecol])
            molsmi.append(x)
        i += 1
    # MORGAN Fingerprints.
    fps = [AllChem.GetMorganFingerprint(x,2,useFeatures=TURE,nBits=1024) for x in molsmi]
    # Note above: multi parameters can be used to generate E/FCFP.
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) # how to name the col?
    df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    df['CHEMBL'] = df.index
    return(df)

def get_fingerprint_from_DataFrame(chem_smile,fpfunc):
    molsmitmp = [Chem.MolFromSmiles(x) for x in chem_smile['smiles']]
    i = 0
    molsmi = []
    for x in molsmitmp:
        if x is not None:
            x.SetProp("_Name",chem_smile['compound'][i])
            molsmi.append(x)
        i += 1
    fps = [fpfunc(x) for x in molsmi]
    # Note above: multi parameters can be used to generate E/FCFP.
    fpsmat = np.matrix(fps)
    df = DataFrame(fpsmat,index = [x.GetProp("_Name") for x in molsmi]) 
    df.insert(0,'chembl',df.index)
    df.insert(1,'smiles',[Chem.MolToSmiles(x) for x in molsmi])
    #df['SMILES'] = [Chem.MolToSmiles(x) for x in molsmi]
    #df['CHEMBL'] = df.index
    return(df)

def main(inputfile,namecol, smicol,outputfile):
    phychemFileStream(inputfile,namecol,smicol,outputfile,chunknum=1000)
#   Some compounds might lose because of the illegal format of SMILES.
#   Thus we need more careful about the results and retrieve those lost compounds.


if __name__ == "__main__":
    inputfile = "chemble2cid2smiles.txt"
    outputfile = "phychemRDKit_smiles2chembl.txt"
    namecol = 1
    smicol = 2
    main(inputfile,namecol,smicol,outputfile)
    


    
