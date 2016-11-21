# retrieve uniprot information.
# szu
# 2016-11-21

from bioservices import UniProt
from matplotlib import pyplot as plt


def quick_getprotinfo(protlist):
    """get protein information from uniprot database
    uniprot(http://www.uniprot.org) based on the package of
    bioservices.

    input ::=
    protlist: list of proteins
    idnm: type of protein names, such as AC.
    output::
    dict of protein information::
    Entry name; Gene names; Length; Organism; Protein names; Status.
    """
    u = UniProt(verbose=False)
    return u.quick_search(protlist)


# if "__name__" == "__main__":
u = UniProt(verbose=False)
data = u.search(
    "zap70+taxonomy:9606",
    frmt="tab",
    limit=3,
    columns="entry name, length, id, genes")
print data

res = u.search(
    "DNMT1_HUMAN",
    frmt="tab",
    columns="entry name, protein names, pathway, comments")
print(res)

res = u.quick_search("DNMT1_HUMAN")
getreskey = res.keys()
res[getreskey]['Gene names']

df = u.get_df("GALK1_HUMAN")
df['Length'].hist()
plt.show()
