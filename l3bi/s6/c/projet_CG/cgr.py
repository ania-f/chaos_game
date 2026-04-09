"""
cgr.py
Calcul des coordonnées CGR et de la matrice de signature génomique.
"""

import numpy as np

# Coordonnées des 4 coins du carré CGR
CORNERS = {
    'A': (0, 1),
    'T': (1, 1),
    'G': (1, 0),
    'C': (0, 0)
}

def generate_cgr_coordinates(sequence):
    """
    Génère les coordonnées CGR pour chaque nucléotide de la séquence.

    Paramètres
    ----------
    sequence : str
        Séquence ADN.

    Retourne
    --------
    coordinates : list of tuple (x, y)
        Une coordonnée par nucléotide.
    """
    x, y = 0.5, 0.5
    coordinates = []
    for nuc in sequence:
        if nuc in CORNERS:
            cx, cy = CORNERS[nuc]
            x = (x + cx) / 2
            y = (y + cy) / 2
            coordinates.append((x, y))
    return coordinates


def generate_signature_matrix(sequence, k):
    """
    Génère la matrice de fréquences CGR pour des k-mers (k de 1 à 4).
    La matrice est de taille 2^k x 2^k, normalisée par le nombre total de k-mers.

    Paramètres
    ----------
    sequence : str
        Séquence ADN.
    k : int
        Taille des k-mers (1 à 4).

    Retourne
    --------
    counts : np.ndarray
        Matrice de fréquences normalisées de taille (2^k, 2^k).
    """
    if k < 1 or k > 4:
        raise ValueError("k doit être compris entre 1 et 4.")

    n_bins = 2 ** k
    counts = np.zeros((n_bins, n_bins), dtype=np.float64)

    coordinates = generate_cgr_coordinates(sequence)

    # La coordonnée à l'index (k-1+i) correspond au k-mer finissant à la position i
    for idx in range(k - 1, len(coordinates)):
        x, y = coordinates[idx]
        i = min(int(x * n_bins), n_bins - 1)
        j = min(int(y * n_bins), n_bins - 1)
        counts[j, i] += 1

    # Normalisation
    total = counts.sum()
    if total > 0:
        counts /= total

    return counts