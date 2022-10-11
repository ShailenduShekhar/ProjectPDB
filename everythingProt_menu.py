#!/usr/bin/python3.8

import sys
from everythingProt_utilities import CleanPDB, ReadPDB


if len(sys.argv) == 1:
    print("Use the '-h' option to open the help page.")

#shows the help page
elif len(sys.argv) == 2 and sys.argv[1] == "-h":
    print("\n\n***  This program generates important information about a given pdb file (which is still in the making) ***")
    print("""


                -c: Displays the number of chains a pdb file has, along with the length of each chain.
                    USAGE: eprot -c input_filename 

                    (While using this option, please make sure that the pdb file has been cleaned, as this tells the number of 
                        C-alpha atoms by reading the ATOM part of the PDB file, and not the SEQRES part; To clean your PDB file
                        please use the option given below)

            -clean: removes ANISOU, Hydrogen atoms, HETATMs (not part of chains) and Alternate Atoms from the pdb file
                    USAGE: eprot -clean input_filename output_filename

                    SUB-OPTIONS:

                        -clean an : Removes only the ANISOU entries
                        -clean hy : Removes only the entries containing Hydrogen atoms
                        -clean het: Removes only the entries containing HETATMs (excluding those that make up the chains)
                        -clean alt: Removes only the ATOM entries that occupy less than the maximum occupancy found among all 
                                    the alternate forms

              -nsa: detect or mutate the non-standard amino acids present in a given PDB file
                    USAGE:
                        
                        eprot -nsa input-file-name             (prints the name and frequency of non-standards amino acids, if any)
                        eprot -nsa input-file-name output-file (mutates the non-standard amino acid(s) to the nearest standard 
                                                                amino acid, if any)

                                      <=====  THANK YOU   =====> 


                                      """)
#shows information about the number of chains and their lengths
elif len(sys.argv) == 3 and sys.argv[1] == '-c':
    ob = ReadPDB(sys.argv[2])
    ob.chainsInfo()

elif len(sys.argv) == 3 and sys.argv[1] == "-nsa":
    ob = ReadPDB(sys.argv[2])
    ob.detectNSA()

elif len(sys.argv) == 4 and sys.argv[1] == "-clean":
    ob = CleanPDB(sys.argv[2],sys.argv[3])

elif len(sys.argv) == 5 and sys.argv[1] == "-clean" and sys.argv[2] == "an":
    ob = CleanPDB()
    ob.cleanANISOU(sys.argv[3],sys.argv[4])

elif len(sys.argv) == 5 and sys.argv[1] == "-clean" and sys.argv[2] == "hy":
    ob = CleanPDB()
    ob.cleanH(sys.argv[3],sys.argv[4])

elif len(sys.argv) == 5 and sys.argv[1] == "-clean" and sys.argv[2] == "het":
    ob = CleanPDB()
    ob.cleanHETATM(sys.argv[3],sys.argv[4])

elif len(sys.argv) == 5 and sys.argv[1] == "-clean" and sys.argv[2] == "alt":
    ob = CleanPDB()
    ob.cleanAlt(sys.argv[3],sys.argv[4])

else: print("\nERROR: No such option/utility\n")
