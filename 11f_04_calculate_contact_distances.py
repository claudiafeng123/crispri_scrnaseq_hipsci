#!/usr/bin/env python3
"""
Calculate the contact distances between chains in a PDB and write in TSV format
"""
import argparse
import sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def main():
    """
    Main
    """
    args = parse_args()

    pdb_id = Path(args.pdb_path).stem

    parser = PDBParser()
    pdb = parser.get_structure(pdb_id, args.pdb_path)[0] # Assume one model

    print("chain1", "residue1", "position1", "chain2", "residue2", "position2", "distance",
          sep="\t", file=sys.stdout)
    for chain1 in pdb:
        for chain2 in pdb:
            for residue1 in chain1:
                resname1 = seq1(residue1.resname)
                if resname1 == "" or not residue1.id[0] == " ":
                    continue

                for residue2 in chain2:
                    resname2 = seq1(residue2.resname)
                    if resname2 == "" or not residue2.id[0] == " ":
                        continue

                    try:
                        print(chain1.id, resname1, residue1.id[1],
                              chain2.id, resname2, residue2.id[1],
                              np.linalg.norm(residue1["CA"].coord - residue2["CA"].coord),
                              sep="\t", file=sys.stdout)
                    except KeyError as err:
                        print("KeyError: ", str(err), " for residue ", residue1, " and ", residue2,
                              sep="", file=sys.stderr)

def arg_parser():
    """
    Construct argument parser
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("pdb_path", metavar="P", help="PDB file path")

    return parser

def parse_args():
    """
    Parse and validate script arguments
    """
    args = arg_parser().parse_args()

    return args

if __name__ == "__main__":
    main()