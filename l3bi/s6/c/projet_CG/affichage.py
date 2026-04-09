"""
affichage.py
Visualisation de la matrice de signature génomique avec matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_signature_matrix(counts, gbk_file, k):
    """
    Affiche la matrice de signature génomique sous forme d'image.

    Paramètres
    ----------
    counts : np.ndarray
        Matrice de fréquences normalisées (2^k x 2^k).
    gbk_file : str
        Nom du fichier GenBank (pour le titre).
    k : int
        Taille des k-mers.
    """
    n_bins = 2 ** k
    plt.figure(figsize=(8, 7))
    cmap = plt.get_cmap('gist_yarg')
    plt.imshow(counts, cmap=cmap, origin='lower')

    # Affichage des valeurs dans chaque case
    fontsize = max(5, 10 - k)
    for i in range(n_bins):
        for j in range(n_bins):
            valeur = counts[i, j]
            text_color = 'black' if valeur < np.max(counts) / 2 else 'white'
            plt.text(j, i, f"{valeur:.4f}", ha='center', va='center',
                     color=text_color, fontsize=fontsize)

    plt.colorbar(label='Fréquence relative', orientation='vertical')
    plt.title(f"Signature génomique avec {k}-mers ({gbk_file})")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.show()


def save_matrix(counts, k, gbk_file, output_file):
    """
    Sauvegarde la matrice de signature dans un fichier texte.

    Paramètres
    ----------
    counts : np.ndarray
        Matrice de fréquences normalisées.
    k : int
        Taille des k-mers.
    gbk_file : str
        Nom du fichier GenBank source.
    output_file : str
        Chemin du fichier de sortie.
    """
    import os
    if os.path.exists(output_file):
        print("Attention: fichier existant, les données seront ajoutées.")

    with open(output_file, 'a') as f:
        f.write("\n----------------------------------------\n")
        f.write(f"Signature génomique avec des {k}-mers (fichier : {gbk_file})\n")
        np.savetxt(f, counts, fmt='%.6f')
    
    print(f"Matrice sauvegardée dans {output_file}")