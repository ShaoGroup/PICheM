# retrieve uniprot information.
# szu
# 201

import uniprot
import pprint

if "__name__" == "__main__":
    seqids = "NP_000508.1 NP_001018081.3".split()
    pairs = uniprot.batch_uniProt_id_mapping_pairs('P_REFSEQ_AC', 'ACC',
                                                   seqids)
    pprint.pprint(pairs, indent=2)
    uniprot_seqids = 'A0QSU3 D9QCH6'.split()
    uniprot_data = uniprot.batch_uniprot_metadata(uniprot_seqids, 'cache')
    pprint.pprint(mapping, indent=2)
