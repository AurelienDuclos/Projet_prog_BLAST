import math

# Création du dictionnaire kmer pour ref et requete
requete = "TTTTTTTTTMKWVTFISLLLARTMKERAE"
database = [
    "AAAAAAAAAMKWVTFISGGGARTMKERAE",
    "MKWVTFISLLLARTMKERAE",
    "GGGGGGGGGMKWVTFISSSSARTMKERAE",
    "TTTTTAEREATDGARAKMGGMAAKED"
]

# Taille des k-mers
k = 4

""" 
Fonction qui coupe une séquence en k-mers.
Prends en entrée une séquence et retourne un dictionnaire,
contenant les k-mers et leurs positions sur la séquence.
"""
def parsing(seq):
    kmers = {}
    for aa in range(len(seq)-k+1):
        kmer = seq[aa:aa+k]
        if kmer not in kmers:
            kmers[kmer] = []
        kmers[kmer].append(aa)
    return kmers


"""
Construction de la matrice blosum62
"""
blosum62 = {}
aa = "ARNDCQEGHILKMFPSTWYV"

matrice = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1,-1,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4, 0,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2,-1,-1, 0,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4]
]

for i, aa1 in enumerate(aa):   # i correspond à la position de l'aa
    blosum62[aa1] = {}
    for j, aa2 in enumerate(aa):
        blosum62[aa1][aa2] = matrice[i][j]


"""
Fonction qui calcule le score de chaque aa.
Prends en entrée les kmers des deux séquences et retourne le score des kmers
"""
def calcul_score(kmer1, kmer2):
    score = 0
    for aa1, aa2 in zip(kmer1, kmer2):
        score += blosum62[aa1][aa2]
    return score



"""
Fonction qui donne une liste de mots similaire au kmer
prends en entrée un kmer et un seuil T et renvoie une liste de mots
"""
def mots_voisins(kmer, T):
    voisins = []
    for i in range(len(kmer)):
        for aa in blosum62:
            mot = kmer[:i] + aa + kmer[i+1:]
            score = calcul_score(kmer, mot)
            if score >= T:
                voisins.append(mot)

    return voisins



"""
Fonction qui trouve les hits (k-mers > seuil T) entre la séquence requete et referente.
Prends en entrée les dictionnaires de kmers, ainsi qu'un seuil T et retourne une liste contenant
les hits et leurs positions sur les deux séquences.
"""
def hits(dico_req, dico_ref, T):
    hits = []
    for kmer in dico_req:
        voisins = mots_voisins(kmer, T)

        for mot in voisins:
            if mot in dico_ref:

                for pos_req in dico_req[kmer]:
                    for pos_ref in dico_ref[mot]:
                        hits.append((kmer, mot, pos_req, pos_ref))

    return hits


"""
Fonction qui étend les hits.
Prends en entrée un hit et le X-dropoff et retourne
un tuple contenant les hits étendues, leurs positions sur les
deux séquences et leur score.
"""
def extension(hit, X, requete, reference):
    kmer_req, kmer_ref, pos_req, pos_ref = hit
    
    kmer_req = kmer_req
    kmer_ref = kmer_ref

    seed_score = calcul_score(kmer_req, kmer_ref)

    best_score_droite = seed_score
    best_score_gauche = seed_score

    debut_req = pos_req
    debut_ref = pos_ref

    fin_req = pos_req + k - 1
    fin_ref = pos_ref + k - 1

    # Extension à gauche
    score = seed_score
    i = pos_req - 1
    j = pos_ref - 1

    while i >= 0 and j >= 0:
        score += blosum62[requete[i]][reference[j]]
        if score > best_score_gauche:
            best_score_gauche = score

            debut_req = i
            debut_ref = j
            
        if best_score_gauche - score > X:
            break
        i -= 1
        j-= 1

    # Extension à droite
    score = seed_score  # Réinitialise le score
    i = pos_req + k
    j = pos_ref + k

    while i < len(requete) and j < len(reference):
        score += blosum62[requete[i]][reference[j]]
        if score > best_score_droite:
            best_score_droite = score

            fin_req = i
            fin_ref = j

        if best_score_droite - score > X:
            break
        i += 1
        j+= 1

    kmer_req = requete[debut_req:fin_req + 1]
    kmer_ref = reference[debut_ref:fin_ref + 1]

    # Score HSP
    score = calcul_score(kmer_req, kmer_ref)  # Réinitialise le score
    msp = best_score_gauche + best_score_droite - seed_score
    if msp > 30:
        return kmer_req, kmer_ref, debut_req, fin_req, debut_ref, fin_ref, msp

    return None


"""
Calcule la taille totale de la database
"""
def taille_database(database):
    total = 0
    for sequence in database:
        total += len(sequence)
    return total


"""
Calcule l'E-value à partir du score MSP
"""
def calcul_evalue(msp, m, n, K, lambda_):
    e_value = K * m * n * math.exp(-lambda_ * msp)
    return e_value
    


"""
Main, affiche les alignements, les positions, leur score MSP et leurs statistiques
"""
T = 17
X = 15

dico_requete = parsing(requete)
resultats = []

for num_sequence, reference in enumerate(database):

    print("\n===================================")
    print("Séquence", num_sequence + 1)
    print("===================================")

    dico_reference = parsing(reference)
    liste_hits = hits(dico_requete, dico_reference, T)

    hit_ext = []
    for hit in liste_hits:
        result = extension(hit, X, requete, reference)
        if result is not None:
            if result not in hit_ext:
                hit_ext.append(result)

    if hit_ext:
        best_msp = hit_ext[0]
        for msp in hit_ext:
            if msp[-1] > best_msp[-1]:
                best_msp = msp
        resultats.append((num_sequence + 1, best_msp))

    print(hit_ext)


m = len(requete)
n = taille_database(database)

k = 0.134
lambda_ = 0.318

print("\n========== STATISTIQUES ==========")

for num_sequence, msp in resultats:
    score_msp = msp[-1]
    e_value = calcul_evalue(score_msp, m, n, k, lambda_)
    print("Séquence :", num_sequence)
    print("MSP :", score_msp)
    print("E-value :", e_value, "\n")
 