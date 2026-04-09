# Projet Chaos Game - Représentation CGR de génomes

## Auteure
Ania FRIK

## Description
Ce projet implémente la technique du Chaos Game Representation (CGR) 
pour analyser et comparer des génomes.

## Fichiers

- `chaos_game.py` : Représentation CGR simple d'un génome
- `sequence.py` : Lecture des fichiers GenBank
- `cgr.py` : Calcul des coordonnées et matrices CGR
- `affichage.py` : Visualisation et sauvegarde
- `genomic_sign.py` : Signature génomique (k-mers de 1 à 4)
- `clustering.py` : Classification hiérarchique et dendrogramme

## Utilisation

```bash
# Partie 1 : CGR simple
python3 chaos_game.py ./genomes/fichier.gbk

# Partie 2 : Signature génomique
python3 genomic_sign.py ./genomes/fichier.gbk <k> <fichier_sortie>

# Partie 3 : Clustering
python3 clustering.py ./genomes/ <k>
```

## Exemple
```bash
python3 genomic_sign.py ./genomes/EU810403.gbk 3 resultats.txt
python3 clustering.py ./genomes/ 5
```

## Référenc
1. Deschavanne et al., Molecular Biology and Evolution, 16, 1391–1399 (1999)
2. Löchel et al., Computational and Structural Biotechnology Journal, 19, 6263-6271 (2021)
