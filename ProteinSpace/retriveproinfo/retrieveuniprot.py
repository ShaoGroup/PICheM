# retrieve uniprot information.
# szu
# 201

from bioservices import UniProt
from matplotlib import pyplot as plt

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

df = u.get_df("GALK1_HUMAN")
df['Length'].hist()
plt.show()
