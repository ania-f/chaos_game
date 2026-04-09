"""
Nom du fichier : clustering.py
Description : Ce fichier contient le code source de l'implémentation d'un clustering hiérarchique basé sur la signature génomique d'une séquence. Il génère le dendrogramme correspondant.
Auteur : Assa DIABIRA & Inès MANOUR
Dernière modification : 22/04/2024
"""
"""
#-----------------------------------------------------------------------------------------------------------------
#                                   Partie 3 :  Classification hiérarchique
#-----------------------------------------------------------------------------------------------------------------

import numpy as np
import sys
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import os

def read_sequence_from_gbk(gbk_file):
    """
    '''Lecture du fichier GenBank et extraction de la séquence d'ADN.'''
    """
    sequence = ""
    with open(gbk_file, "r") as file:
        for line in file:
            if line.startswith("ORIGIN"):
                break
        for line in file:
            if line.startswith("//"):
                break
            filtered_line = ''.join(filter(str.isalpha, line.strip().upper()))
            sequence += ''.join(filter(lambda base: base in ['A', 'C', 'G', 'T'], filtered_line))
    return sequence

def generate_signature_matrix(sequence, k):
    """
    '''Génère une matrice de signature génomique à partir de la séquence d'ADN.'''
    """
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    signature_matrix = np.zeros((4, 4), dtype=np.int64)
    bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for kmer in kmers:
        for i in range(len(kmer) - 1):
            current_base = kmer[i]
            next_base = kmer[i + 1]
            signature_matrix[bases[current_base], bases[next_base]] += 1
    return signature_matrix

def calculate_similarity_matrix(signature_matrices):
    """
   ''' Calcule une matrice de similarité entre les séquences d'ADN à partir de leurs signatures génomiques.'''
    """
    n = len(signature_matrices)
    similarity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            similarity_matrix[i, j] = calculate_similarity(signature_matrices[i], signature_matrices[j])
    return similarity_matrix

def calculate_similarity(matrix1, matrix2):
    """
    '''Calcule la similarité entre deux matrices de signature génomique.''''
    """
    similarity = np.corrcoef(matrix1.flatten(), matrix2.flatten())[0, 1]
    return similarity

def hierarchical_clustering(similarity_matrix, gbk_files):
    """
    '''Effectue le clustering hiérarchique à partir de la matrice de similarité et affiche le dendrogramme.'''
    """
    linkage_matrix = linkage(similarity_matrix, method='average')
    plt.figure(figsize=(10, 8))

    # print("Dimensions de la matrice de liaison:", linkage_matrix.shape) vérifications

    dendrogram(linkage_matrix, labels=gbk_files, leaf_rotation=90)
    plt.title('Dendrogramme de la Classification hiérachique')
    plt.xlabel('Séquences GenBank testées')
    plt.ylabel('Similarité')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    directory = sys.argv[1]
    k = 5
    signature_matrices = []
    gbk_files = []
    for file in os.listdir(directory):
        if file.endswith(".gbk"):
            gbk_files.append(file)
            sequence = read_sequence_from_gbk(os.path.join(directory, file))
            if sequence:
                signature_matrix = generate_signature_matrix(sequence, k)
                signature_matrices.append(signature_matrix)
            else:
                print("Impossible de générer une séquence à partir du fichier:", file)


    print("Nombre de fichiers GenBank trouvés:", len(gbk_files))
    print("Nombre de matrices de signature génomique générées:", len(signature_matrices))

    similarity_matrix = calculate_similarity_matrix(signature_matrices)
    hierarchical_clustering(similarity_matrix, gbk_files)"""

    """ version 2 """

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

def read_sequence_from_gbk(gbk_file):
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

def generate_cgr_coordinates(sequence):
    corners = {'A': (0, 1), 'T': (1, 1), 'G': (1, 0), 'C': (0, 0)}
    x, y = 0.5, 0.5
    coordinates = []
    for nuc in sequence:
        if nuc in corners:
            cx, cy = corners[nuc]
            x = (x + cx) / 2
            y = (y + cy) / 2
            coordinates.append((x, y))
    return coordinates

def generate_signature_matrix(sequence, k):
    n_bins = 2 ** k
    counts = np.zeros((n_bins, n_bins), dtype=np.float64)
    coordinates = generate_cgr_coordinates(sequence)
    for idx in range(k - 1, len(coordinates)):
        x, y = coordinates[idx]
        i = min(int(x * n_bins), n_bins - 1)
        j = min(int(y * n_bins), n_bins - 1)
        counts[j, i] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts

def calculate_distance_matrix(signature_matrices):
    n = len(signature_matrices)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            v1 = signature_matrices[i].flatten()
            v2 = signature_matrices[j].flatten()
            dist = np.sqrt(np.sum((v1 - v2) ** 2))
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    return distance_matrix

def hierarchical_clustering(distance_matrix, labels):
    condensed = squareform(distance_matrix)
    linkage_matrix = linkage(condensed, method='average')
    plt.figure(figsize=(12, 7))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=45, leaf_font_size=10)
    plt.title('Dendrogramme - Classification hiérarchique (CGR k-mers)', fontsize=14)
    plt.xlabel('Génomes', fontsize=12)
    plt.ylabel('Distance CGR', fontsize=12)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python clustering.py <dossier_genbank> [k]")
        sys.exit(1)

    directory = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) >= 3 else 5

    signature_matrices = []
    gbk_files = []

    for file in sorted(os.listdir(directory)):
        if file.endswith(".gbk"):
            path = os.path.join(directory, file)
            sequence = read_sequence_from_gbk(path)
            if sequence:
                print(f"  ✓ {file} ({len(sequence)} pb)")
                sig = generate_signature_matrix(sequence, k)
                signature_matrices.append(sig)
                gbk_files.append(file)
            else:
                print(f"  ✗ {file} : séquence vide")

    print(f"\n{len(gbk_files)} génomes chargés, k={k}")

    
    
    distance_matrix = calculate_distance_matrix(signature_matrices)
    hierarchical_clustering(distance_matrix, gbk_files)