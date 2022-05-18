# -*- coding: utf-8 -*-
from Bio import SeqIO  # Tip: This module might be useful for parsing...


############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    PARSER = SeqIO

    def __init__(self, path, frmt='uniprot-xml'):
        """
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        """
        self.record = self.PARSER.read(path, frmt)  # Parse the XML file once and re-use it in the functions below

    # 2.2 SwissProt Identification
    def get_sp_identification(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                identification: tuple consisting of in this order
                                    1. string: Unique SwissProt identifier for the given xml file
                                    2. string: Primary protein name
                                    3. string: Primary gene name
        """
        return (
            self.record.id,
            self.record.annotations['recommendedName_fullName'][0],
            self.record.annotations['gene_name_primary'],
        )

    # 2.3 SwissProt Sequence Information
    def get_sp_sequence_info(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                information: tuple consisting of in this order
                                    1. str: sequence of the UniProt entry
                                    2. int: sequence length of the UniProt entry
                                    3. int: sequence mass of the UniProt entry
        """
        return (
            str(self.record.seq),
            self.record.annotations['sequence_length'],
            self.record.annotations['sequence_mass'],
        )

    # 2.4 Organism 
    def get_organism(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        """
        return self.record.annotations['organism']

    # 2.5 Localizations
    def get_localization(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the subcellular localization as stcated in the
                corresponding field.
                Return value has to be a list of strings.
        """
        return self.record.annotations['comment_subcellularlocation_location']

    # 2.6 Cross-references to PDB
    def get_pdb_support(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        """
        prefix = 'PDB:'
        return [r.removeprefix(prefix) for r in self.record.dbxrefs if r.startswith(prefix)]


def main():
    print('SwissProt XML Parser class')
    return None


if __name__ == '__main__':
    main()
