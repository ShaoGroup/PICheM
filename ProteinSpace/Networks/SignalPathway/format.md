##**Egde table**
**node1**

**node2**

**Subtype1**	_activation;binding/association;dissociation;indirect effect;NULL_

**Subtype2**	_dephsophorylation;phosphorylation;NULL_

**EdgeType**	_-->;---;..>;--|;-+-;NULL_

**Source**	_hsa04110;CST00_


##**Node table**
**node**

**type**	_protein;DNA/RNA;compound;pathway_

**gene**(for protein node only)

**cellular_component**	_cytosol;nucleus;plasma membrane;mitochondrial membrane;ER membrane;else?_


##**Source table**
**source**	_KEGG;CST_

**PathwayID**(correspond to Source in edge table)

**PathwayName**

**Label**

**Link**

**comment**


#_NOTE_
1. All tables shall be comma seperated files (.csv);
2. About EdgeType (for KEGG orginated edges only):

>--> activation

>--| inhibition

>..> indirect effect

>--- binding

>-+- association 

