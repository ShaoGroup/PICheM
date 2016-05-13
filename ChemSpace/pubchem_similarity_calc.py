# -*- coding: utf-8 -*- 

# Downdload larrge-scale compounds' fingerprints from the PubChem Database.
# ?Calculate the similarities 
# Songpeng Zu
# 2015-11-10 

# import package
import pubchempy as pcb
import chemfp.bitops as chembit
import argparse # Temp for args.

# Functions
def chunks_get_compounds(l,n):
    """ Yield successive n-sized chunks from l. (From Stack Overflow)"""
    for i in xrange(0,len(l),n):
        try:
            yield pcb.get_compounds(l[i:i+n])
        except Exception, e:
            print e
            pass

def read_hexfile(infile="seed_drug_cid2fingerprint"):
    """ Read all the lines into a list. 

        First column as the cids, Second column as the hex string.
    """
    with open(infile,'r') as fr:
        cids = [x.strip().split()[0] for x in fr]
    with open(infile,'r') as fr:
        fps = [x.strip().split()[1] for x in fr]
    return cids, fps

def download_PubChemFinger(line_start_num = 0,line_end_num = 1000,block_num=100,
                           inputfile="zinc2cid_lead_uniq",
                           outfilename="zinc_cid2pubchemfp",colnum=2,cid_col_index=1):
    """Directly download the fingerprints from PubChem.
    
    Note: if one column per line (colnum==1), we treat it as CID directly;
          if multiple column per line (colnum >1), we use the cid_col_index for the column specific for CID(PubChem ID).
          some of them return more than one fingerprints, and we didnot consider this situation.
    """
    # Get compounds_list
    compounds_list = []
    with open(inputfile,'r') as fr:
        if colnum < 2:
            compounds_list = [int(line.strip()) for line in fr]
        else:
            compounds_list = [int(line.strip().split()[cid_col_index]) for line in fr]
    # Write the results.
    with open(outfilename+"_"+str(line_start_num)+"_"+str(line_end_num),'w') as fw:
        compounds_list = compounds_list[line_start_num:line_end_num] # Resize compounds_list
        for compounds in chunks_get_compounds(compounds_list,block_num):
            try:
                tmp_write_list = ['\t'.join([str(compound.cid),compound.fingerprint]) for compound in compounds]
                fw.write('\n'.join(tmp_write_list))
            except Exception,e: # Ignore the possible errors.
                print e # print the error, then continue
                continue

def tanimoto_by_MattSwain(compound1,compound2):
    """This function is provided from pubchempy by MattSwain.
       
       It's not fast. Using the hex_tanimoto form chemfp.bitops instead.
    """
    fp1 = int(compound1.fingerprint, 16) # as HEX
    fp2 = int(compound2.fingerprint, 16)  
    fp1_count = bin(fp1).count('1') # binary
    fp2_count = bin(fp1).count('1')
    both_count = bin(fp1 & fp2).count('1')
    return float(both_count) / (fp1_count + fp2_count - both_count)

def cal_PubChem_Similarity_one2all(one_zinc_fp_hex,fp_hex_drug):
    """Calculate the similarites for one zinc compounds against all compounds.
    """
    tmp_array = [round(chembit.hex_tanimoto(one_zinc_fp_hex,fp_drug),3)
            for fp_drug in fp_hex_drug]
    return '\t'.join([str(x) for x in tmp_array])
    
def calc_PubChem_Similarity(infilenm_zinc,infilenm_drug,outfilenm):
    """Calculate the PubChem Similarites between compounds from zinc and drugbank.
    """
    cid_zinc, fp_hex_zinc = read_hexfile(infilenm_zinc)
    cid_drug, fp_hex_drug = read_hexfile(infilenm_drug)
    with open(outfilenm + '_cid','w') as fw_cid:
        fw_cid.write('\n'.join([str(cid) for cid in cid_zinc]))
    with open(outfilenm,'w') as fw:
        fw.write('\n'.join(
            [cal_PubChem_Similarity_one2all(fp_zinc,fp_hex_drug)
             for fp_zinc in fp_hex_zinc]))

def download_fp():
    parser = argparse.ArgumentParser()
    parser.add_argument("line_start",type=int,help="The start line for reading")
    parser.add_argument("line_end",type=int,help="The end line for reading")
    args = parser.parse_args()
    download_PubChemFinger(args.line_start,args.line_end) # Set the line start and line end.
    
def get_similarity_matrix():
    parser = argparse.ArgumentParser()
    parser.add_argument("zinc_file_name",type=str,help="The input file name for zinc")
    parser.add_argument("drugbank_file_name",type=str,help="The input file name for drugbank")
    parser.add_argument("output_file",type=str,help="The output file name for the similarity matrix")
    args = parser.parse_args()
    calc_PubChem_Similarity(args.zinc_file_name,args.drugbank_file_name,args.output_file)
    

if __name__ == "__main__":
#    download_fp()
    get_similarity_matrix()
