"""
genomic_sign.py
Programme principal - Signature génomique par CGR.
Utilisation : python3 genomic_sign.py <fichier.gbk> <k> <fichier_sortie>
"""

import sys
from sequence import read_sequence_from_gbk
from cgr import generate_signature_matrix
from affichage import plot_signature_matrix, save_matrix


if __name__ == "__main__":

    # Vérification des arguments
    if len(sys.argv) != 4:
        print("Erreur: Nombre incorrect d'arguments.")
        print("Usage: python3 genomic_sign.py <fichier_genbank> <k> <fichier_sortie>")
        sys.exit(1)

    gbk_file = sys.argv[1]
    output_file = sys.argv[3]

    # Vérification de k
    try:
        k = int(sys.argv[2])
        if k < 1 or k > 4:
            raise ValueError
    except ValueError:
        print("Erreur: k doit être un entier entre 1 et 4.")
        sys.exit(1)

    # Lecture de la séquence
    print(f"Lecture du fichier GenBank : {gbk_file}")
    sequence = read_sequence_from_gbk(gbk_file)

    if not sequence:
        print("Erreur: séquence vide ou fichier inaccessible.")
        sys.exit(1)

    print(f"Séquence lue : {len(sequence)} nucléotides")

    # Calcul de la matrice CGR
    counts = generate_signature_matrix(sequence, k)

    # Affichage
    plot_signature_matrix(counts, gbk_file, k)

    # Sauvegarde
    save_matrix(counts, k, gbk_file, output_file)