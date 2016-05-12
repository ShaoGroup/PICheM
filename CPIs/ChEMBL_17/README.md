## ChEMBL Data Version 17
#### Collected by: Songpeng Zu, Jianzhou Zhou
#### Last update time: May, 2015

Here we collected Compound-Protein Interactions (CPIs) from the ChEMBL database,
version 17. The principles for these data are:
1. Targets labeled as "SINGLE PROTEIN"  
2. Assays labeled as "B" (Binding).  
3. Confidence Scores no less than 7.  
4. IC50/EC50/Ki/Kd treated no distinctions, and the unit is nM.  
See the file "compounds_targets_interaction_7.txt" for details.


Proteins also have the organisms ("organism.txt") and classification
("Protein2Class_Chembl.txt") information. 

Here we upload a R script ("CPIs_Data.R") as a demo to deal with the data
set. And the human proteins-compound interactions with protein classification
information are uploaded as "CPIsByR.txt". 
