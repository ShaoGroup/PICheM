#!/usr/bin/python

# Extract drugs and experimental compounds.
# Songpeng Zu
# 2014-12-30

import os
import xml.sax

class DrugBankHandler(xml.sax.ContentHandler):
    def __int__(self):
        self.CurrentData = ""
        self.DrugBankID = ""
        self.smallmolecule = ""
        self.drugbank = ""
        self.parentContent = ""
        self.SMILES = ""
    def startElement(self,tag,attributes):
        self.CurrentData = tag
        if tag == "drug" and  "type" in attributes:
            if attributes["type"] == "small molecule":
                self.smallmolecule = "small molecule"
        if tag == "drugbank-id" and 'primary' in attributes:
            self.drugbank = True
    def endElement(self,tag):
        if self.CurrentData == "drugbank-id":
            if self.drugbank and "DB0" in self.DrugBankID:
#                with open("DrugBankID2SMILES.txt","a") as fid:
#                    fid.write(self.DrugBankID+"\t")
                print "DrugBank-ID:",self.DrugBankID+"\n"
        elif self.CurrentData == "value" and self.parentContent == "SMILES":
            if self.smallmolecule and self.drugbank:
                with open("DrugBankID2SMILES.txt","a") as fid:
                    fid.write(self.DrugBankID+"\t"+self.SMILES+"\n")
                self.smallmolecule = ""
                self.drugbank = ""
        self.CurrentData = ""
    def characters(self,content):
        if self.CurrentData == "kind":
            self.parentContent = content
        elif self.CurrentData == "value":
            self.SMILES = content
        elif self.CurrentData == "drugbank-id" and "DB0" in content:
            self.DrugBankID = content

if __name__ == "__main__":
    try:
        os.remove("DrugBankID2SMILES.txt")
    except OSError:
        print "No file has been written."
    finally:
        print "Clean up the written file."

    parser = xml.sax.make_parser()
    parser.setFeature(xml.sax.handler.feature_namespaces,0)
    Handler = DrugBankHandler()
    parser.setContentHandler(Handler)
    parser.parse("drugbank.xml")

    print "The job has been done."



