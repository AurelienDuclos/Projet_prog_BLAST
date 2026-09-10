# Projet BLAST en python

## Fonctionnalités
- Découpage des séquences en k-mers
- Génération des mots voisins avec BLOSUM62
- Recherche des hits
- Extension des hits avec X-dropoff
- Calcul du score MSP
- Calcul de l'E-value

## Exécution
Si vous souhaitez avoir les résultats dans le terminal:
python3 DUCLOS_code.py

Si vous souhaitez avoir les résultats dans un fichier:
python3 DUCLOS_code.py > resultats.txt

## Résultats
Vous obtenez d'abord les alignements pour chaque séquence avec leurs positions et le score MSP
Puis vous avez la E-value du meilleur MSP entre les différents alignements de chaque séquence.
