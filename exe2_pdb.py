# -*- coding: utf-8 -*-
from Bio.Data.SCOPData import protein_letters_3to1
from Bio.PDB.MMCIFParser import MMCIFParser  # Tip: This module might be useful for parsing...
import numpy as np


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    CIF_PARSER = MMCIFParser()  # parser object for reading in structure in CIF format

    def __init__(self, path):
        """
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...
        """
        self.structure = MMCIFParser().get_structure('7AHL', path)

    def get_chain(self, chain_id):
        return next(c for c in self.structure.get_chains() if c.id == chain_id)

    # 2.8 Chains    
    def get_number_of_chains(self):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        """
        return len(list(self.structure.get_chains()))

    # 2.9 Sequence  
    def get_sequence(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        """
        return ''.join([protein_letters_3to1[r.get_resname()] for r in self.get_chain(chain_id).get_residues()
                        if r.get_resname() in protein_letters_3to1])

    # 2.10 Water molecules
    def get_number_of_water_molecules(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        """
        return len([r for r in self.get_chain(chain_id).get_residues() if r.get_resname() == 'HOH'])

    # 2.11 C-alpha distance
    def get_ca_distance(self, chain_id_1, index_1, chain_id_2, index_2):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure

                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.

            Return:
                Return the C-alpha (!) distance between the two residues, described by
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().

            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins
            (Heterodimers) or between different copies of the same protein (Homodimers).
        """
        c1 = self.get_chain(chain_id_1)
        c2 = self.get_chain(chain_id_2)
        return int(1)

    # 2.12 Contact Map  
    def get_contact_map(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet)
                for a given chain in a Biopython.PDB structure as numpy array.
                The values in the matrix describe the c-alpha distance between all residues
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        """
        r = [ri for ri in self.get_chain(chain_id).get_residues() if 'CA' in ri]
        rn = len(r)
        contact_map = np.zeros((rn, rn), dtype=np.float32)
        for i in range(rn):
            for j in range(rn):
                contact_map[i][j] = self.get_ca_distance(chain_id, i + 1, chain_id, j + 1)
        return contact_map.astype(np.int64)  # return rounded (integer) values

    # 2.13 B-Factors    
    def get_bfactors(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue;
                (e.g. the residue was not resolved); insert np.nan for those cases.

                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        """
        r = list(self.get_chain(chain_id).get_residues())
        b = np.array([np.nanmean([a.bfactor for a in ri.get_atoms()]) for ri in r if ri.get_resname() != 'HOH'])
        b -= b.mean()
        b /= b.std()
        return b.astype(np.int64)  # return rounded (integer) values


def main():
    print('PDB parser class.')
    return None


if __name__ == '__main__':
    main()
