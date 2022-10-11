#!/usr/bin/python3.8

import sys
import os
from triplet_codon_dict import aa_dict

class CleanPDB:
    def __init__(self,infile=None,outfile=None):
        if infile != None:
            self.infile = infile
            self.outfile = outfile
            self.cleanEverything()

    def cleanEverything(self):
        self.cleanANISOU(self.infile,"temp1.txt")
        self.cleanH("temp1.txt","temp2.txt")
        self.cleanHETATM("temp2.txt","temp3.txt")
        self.cleanAlt("temp3.txt",self.outfile)
        print(f"\n{self.infile} has been cleaned into {self.outfile}\n")
        os.remove("temp1.txt")
        os.remove("temp2.txt")
        os.remove("temp3.txt")

#removing ANISOU:

    def cleanANISOU(self,infile,outfile):
        print("Removing the ANISOU lines from the file ....................")
        var = open(infile)
        dar = open(outfile,"w")
        line = var.readline()
        while line:
            if "ANISOU" not in line[:6]: dar.write(line)
            line = var.readline()
        var.close()
        dar.close()

#Removing ANISOU works; now removing H-atoms

    def cleanH(self,infile,outfile):
        print("Removing the Hydrogen atoms from the file ..................")
        var = open(infile)
        dar = open(outfile,"w")
        line = var.readline()
        while line:
            if "H" not in line[76:78]: dar.write(line) 
            line = var.readline()
        var.close()
        dar.close()

#Removing H-atoms works; now removing HETATMs from the end of the file

    def cleanHETATM(self,infile,outfile):
        print("Removing the HETATMs from the file .........................")
        var = open(infile)
        dar = open(outfile,"w")
        check = 0
        line = var.readline()
        while line:
            if "TER" in line[:6]: check = 1
            elif check == 1 and "HETATM" not in line[:6]: check = 0
            if check == 1 and line[:6] == "HETATM": break
            dar.write(line)
            line = var.readline()

        dar.write("END")
        dar.close()
        var.close()

#Removing HETATMs from the end of the file is complete; moving onto removing alternate conformations
#13-16 columns of pdb file: Atoms name
#23-26 cols of pdb file: Residue sequence number

    def cleanAlt(self,infile,outfile):
        print("Removing Alternate Atoms from the file .....................")
        def replacing_occu(line):
            li = list(line)
            li[54] = " "
            li[55] = " "
            li[56] = "1"
            li[57] = "."
            li[58] = "0"
            li[59] = "0"
            li[16] = " "
            return "".join(li)

        var = open(infile)
        dar = open(outfile,"w")
        line = var.readline()
        occu_list = []
        prev_name = ""
        line_list = []
        count = 0
        while line:
            if "ATOM" in line[:6] or "HETATM" in line[:6]:
                name = line[12:16] + line[22:26]
                occu = float(line[54:60])
                if occu == 1.00:
                    dar.write(line)
                    line = var.readline()
                    continue
                prev_name = name
                occu_list.append(occu)
                line_list.append(line)
                line2 = var.readline()
                while line2:
                    if "ATOM" in line2[:6] or "HETATM" in line2[:6]:
                        name = line2[12:16] + line2[22:26]
                        occu = float(line2[54:60])
                    else: line = line2; break
                    if name == prev_name:
                        occu_list.append(occu); line_list.append(line2)
                    else:
                        pos = occu_list.index(max(occu_list))
                        to_write = line_list[pos]; to_write = replacing_occu(to_write)
                        dar.write(to_write)
                        occu_list.clear(); line_list.clear()
                        var.seek(var.tell() - len(line2))
                        break
                    line2 = var.readline()
            else: dar.write(line)
            line = var.readline()

class ReadPDB:
    def __init__(self,infile):
        self.infile = infile
        self.preprocessing()

    def preprocessing(self):
        var = open(self.infile)
        self.chain_dict = {}
        self.residue_dict = {}
        self.residue_pos_dict = {}
        self.modified_aa_dict = {}

        prev_res_no = ""
        check = False
        line = var.readline()
        while line:
            if "TER" in line[:6]: 
                check = True
                line2 = var.readline()
                if "HETATM" in line2[:6]: break
                else: line = line2; continue
            if "HETNAM" in line[:6]:
                nsa_id = line[11:14]; chemical_name = line[15:70].strip()
                for residue in list(aa_dict.keys()):
                    if aa_dict[residue] in chemical_name:
                        self.modified_aa_dict[nsa_id] = chemical_name
            if "ATOM" in line[:6] or "HETATM" in line[:6] and "CA" in line[12:16]:
                res_id = line[17:20];
                chainid = line[21]; res_no = line[22:26].strip()
                if res_no != prev_res_no:
                    if self.chain_dict.get(chainid): self.chain_dict[chainid] = self.chain_dict.get(chainid) + [(res_id,res_no)]
                    else: self.chain_dict[chainid] = [(res_id,res_no)]
                    if self.residue_dict.get(res_id): self.residue_dict[res_id] = self.residue_dict.get(res_id) + 1
                    else: self.residue_dict[res_id] = 1
                    #if self.modified_aa_dict.get(res_id): self.modified_aa_dict[res_id] = self.modified_aa_dict.get(res_id) + [res_no]
                    #else: self.modified_aa_dict[res_id] = [res_no]
                    prev_res_no = res_no

            line = var.readline()

    def chainsInfo(self):
        print(f"\nChainId\tChainLen\tStart\tEnd\n")
        for key in list(self.chain_dict.keys()):
            total_len = len(self.chain_dict.get(key))
            print(f"{key}\t{total_len}\t{self.chain_dict.get(key)[0]}\t{self.chain_dict.get(key)[total_len-1]}")

    def detectNSA(self):
        if self.modified_aa_dict:
            location_dict = {}
            for item in list(self.modified_aa_dict.keys()):
                for chain in list(self.chain_dict.keys()):
                    for res,pos in self.chain_dict[chain]:
                        if item == res:
                            if location_dict.get(item): location_dict[item] = location_dict.get(item) + [(chain,pos)]
                            else: location_dict[item] = [(chain,pos)]
            #print(location_dict)
                
            print(f"\nCode\tName\tFreq\tPos\n")
            for item in list(self.modified_aa_dict.keys()):
                print(f"{item}\t{self.modified_aa_dict[item]}\t{self.residue_dict[item]}\t{location_dict[item]}")
            print()
        else:
            location_dict = {}
            for res in list(self.residue_dict.keys()):
                if res not in list(aa_dict.keys()):
                    for chain in list(self.chain_dict.keys()):
                        for resi,pos in self.chain_dict.get(chain):
                            if res == resi: 
                                if location_dict.get(item): location_dict[item] = location_dict.get(item) + [(chain,pos)]
                                else: location_dict[item] = [(chain,pos)]
            if location_dict:
                print(f"\nCode\tFreq\tPos\n")
                for res in list(self.residue_dict.keys()):
                    if res not in list(aa_dict.keys()):
                        print(f"{res}\t{self.residue_dict.get(res)}\t{location_dict.get(res)}")
            else: print("\nNo non-standard amino acid present\n")

            
"""
ob = ReadPDB(sys.argv[1])
#print("modified dictionary")
#print(ob.modified_aa_dict)
#print("chain dictionary")
#print(ob.chain_dict)
#print("residue dictionary")
#print(ob.residue_dict)
ob.detectNSA()
"""
