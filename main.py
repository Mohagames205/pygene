class Strand:
    pass


class DNA(Strand):
    bases = ["A", "C", "T", "G"]

    def __init__(self, dna):
        dna = dna.replace(" ", "")
        for ind_base in list(dna):
            if ind_base not in self.bases:
                raise InvalidStrandException
        self.dna = dna


class RNA(Strand):
    bases = ["A", "C", "T", "U", "G"]
    pretty = ""

    def __init__(self, rna):
        rna = rna.replace(" ", "")
        for ind_base in list(rna):
            if ind_base not in self.bases:
                raise InvalidStrandException

        if "T" in list(rna):
            print("DNA detected, converting to RNA...\n")
            rna = StrandConverter.toRNA(DNA(rna)).rna
        self.rna = rna

    def toProtein(self):

        rna = self.rna

        # look for start codon
        begin = rna.find(StrandConverter.map["START"][0])

        keten = ""
        i = begin + 3
        rnalist = list(rna)
        while i < len(rna):
            if i + 2 < len(rna):
                codon = rnalist[i] + rnalist[i + 1] + rnalist[i + 2]
                self.pretty += "[" + rnalist[i] + rnalist[i + 1] + rnalist[i + 2] + "] "
                if StrandConverter.toAminoAcid(codon) == "STOP":
                    return keten

                keten += StrandConverter.toAminoAcid(codon) + " "
            i += 3
        return keten


class StrandConverter:
    map = {"Phe": ['UUU', 'UUC'],
           "Leu": ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
           "Ile": ["AUU", "AUC", "AUA"],
           "START": ["AUG"],
           "STOP": ["UAA", "UAG", "UGA"],
           "Val": ['GUU', 'GUC', 'GUA', 'GUG'],
           "Ser": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
           "Pro": ['CCU', 'CCC', 'CCA', 'CCG'],
           "Thr": ['ACU', 'ACC', 'ACA', 'ACG'],
           "Ala": ['GCU', 'GCC', 'GCA', 'GCG'],
           "Tyr": ['UAU', 'UAC'],
           "His": ['CAU', 'CAC'],
           "Gln": ['CAA', 'CAG'],
           "Asn": ['AAU', 'AAC'],
           "Lys": ['AAA', 'AAG'],
           "Asp": ['GAU', 'GAC'],
           "Glu": ['GAA', 'GAG'],
           "Cys": ['UGU', 'UGC', 'UGG'],
           "Arg": ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
           "Gly": ['GGU', 'GGC', 'GGA', 'GGG'], }

    @staticmethod
    def toRNA(dna) -> RNA:
        return RNA(dna.dna.replace('T', 'U'))

    @staticmethod
    def toDNA(rna) -> DNA:
        return DNA(rna.rna.replace('U', 'T'))

    @staticmethod
    def toAminoAcid(triplet):
        for codons in StrandConverter.map:
            if triplet in StrandConverter.map[codons]:
                return codons

    @staticmethod
    def detect(strand):
        if "T" in strand:
            return DNA(strand)

        if "U" in strand:
            return RNA(strand)

        return DNA(strand)


class InvalidStrandException(Exception):
    def __init__(self):
        super().__init__('This string is not a valid sequence')


rna = RNA(input("Geef de basensequentie "))
print(rna.toProtein())
