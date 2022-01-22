# Bloom_filter

## Présentation

Ce projet a été réaliser comme épreuve de programmation pour un stage.
Il a pour objectif de réaliser un filtre de Bloom et d'y entrer tous les k-mers dans une séquence ADN donnée (dans un fichier FASTA) puis d'y effectuer des requêtes.

## Compilation

Un Makefile est disponible à la racine du projet. Pour compiler exécuter `make` dans le répertoire du projet.
`make clean` est également disponible pour supprimer l'exécutable produit.

## Exécution 

L'exécution se fait avec une ligne de commande de cette forme : <br> `./hash nom_de_fichier k taille_filtre nombre_hashages nombre_requêtes` <br>
qui est respécifiée en cas de nombre d'arguments incorrect.
Les requêtes et les résultats de ces requêtes sont écrits dans un fichier Requests.txt à la racine du projet en fin de traitement.

### Auteur 

Louis Kurdyk
