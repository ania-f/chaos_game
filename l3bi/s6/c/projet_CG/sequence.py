"""
sequence.py
Lecture et extraction de la séquence ADN depuis un fichier GenBank.
"""

def read_sequence_from_gbk(gbk_file):
    """
    Lit un fichier GenBank et retourne la séquence ADN (uniquement A, C, G, T).
    
    Paramètres
    ----------
    gbk_file : str
        Chemin vers le fichier GenBank.
    
    Retourne
    --------
    sequence : str
        Séquence ADN extraite.
    """
    sequence = ""
    with open(gbk_file, "r") as file:
        in_sequence = False
        for line in file:
            if line.startswith("ORIGIN"):
                in_sequence = True
                continue
            if in_sequence:
                if line.startswith("//"):
                    break
                for char in line.upper():
                    if char in "ACGT":
                        sequence += char
    return sequence