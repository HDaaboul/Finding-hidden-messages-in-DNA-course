"""I put sorted programs from part II bioinformatics course (Stepik):
1-entropy calculator: renders entropy of prob dist of each column in 10*12 motif matrix for an example of Nf-kB consensus recognition site (pretty manual)
2-Motif matrix counter: counts each nt type in each column (assign list of zeros (same len as strings) to each nt as a key in dict, then at each position raise the
count of correct nt at a position of its list same as the position of string reached by for loop)
3-Motif Profile function: a function creates a profile of a string list as a dict with nt in keys and proportion of them in columns; in values
4-Motif Consensus: Renders the most consensus sequence of the motif by taking nt with highes prob in each column in profile
5-Motif Score: gives a score representing the number of mismatched nt among the list of motifs from the consensus (the least the closest to perfect)
6-Pr function: Calculates the probability of the chosen motif to happen according to Profile
7-Profile-Most Probable K-mer: renders the k-mer with highest prob in a string according to a profile.

8-Greedy motif search: applies the greedy algorithm by taking the first k-mer from every string and keep them in 'best motifs' then takes every k-mer from the
first string and compare it with each other string from second to last according to profile of 'motifs' list (with chosen k-mer from first string alone at first)
and when has the profile most probable motif in the second string it appends it to the motifs list and re-do their profile to give pmpp of third one and so on.
each k-mer in the first string is iterated with k-mers in all other strings, and has motifs list with a score. Now, we compare the score of each motifs list to the
default score (score of 'best_motifs' that includes first k-mers. Eventually we find the k-mer in first string that has lowest score, and its motifs list is the
final outcome of the function (as it should be the closest)
"""
####1 entropy calculator

from math import *
def entropy(li):
    totlst = []
    for i in li:
        ent = -(i*log2(i))
        totlst.append(ent)
    totent = sum(totlst)
    return totent
k = 12
li1 = [0.2,0.1,0.7]
li2 = [0.2,0.2,0.6]
li3 = [1]
li4 = [1]
li5 = [0.1,0.9]
li6 = [0.1,0.9]
li7 = [0.1,0.9]
li8 = [0.1,0.4,0.5]
li9 = [0.1,0.1,0.8]
li10 =[0.1,0.2,0.7]
li11 =[0.3,0.3,0.4]
li12 =[0.6,0.4]
limot =[li1,li2,li3,li4,li5,li6,li7,li8,li9,li10,li11,li12]
fin = []
for i in limot:
    fin.append(entropy(i))
print(sum(fin))

####2 Motif matrix counter:

Motifs = []
print("Please enter your lines here and then press Ctrl+D when you finish: ")
while True:
    try:
        line = input()
    except EOFError:
        break
    Motifs.append(line)

def Count(Motifs):
    count = {}
    x = len(Motifs)
    k = len(Motifs[0])
    for sym in "ACTG":
        count[sym] = []
        for r in range(k):
            count[sym].append(0)
    for i in range(x):
        for j in range(k):
            nt = Motifs[i][j]
            count[nt][j] += 1
    return count

####3 Motif Profile: $$needs motif count

def Profile(Motifs):
    profile = Count(Motifs)
    x = len(Motifs)
    k = len(Motifs[0])
    for key in profile:
         for i in range(k):
            profile[key][i] /= x
    return profile
print(Profile(Motifs))

####4 Consensus-rendering func: $$ needs motif count

def Consensus(Motifs):
    cons = ""
    k = len(Motifs[0])
    x = len(Motifs)
    count = Count(Motifs)
    for j in range(k):
        freq = ""
        m = 0
        for nt in "ATGC":
            if count[nt][j]>m:
                m = count[nt][j]
                freq = nt
        cons += freq
    return cons

####5 Score (num of mismatches with consensus: $$needs motif count

def Score(Motifs):
    scr = []
    k = len(Motifs[0])
    x = len(Motifs)
    count = Count(Motifs)
    for j in range(k):
        m = 0
        for nt in "ACTG":
            if count[nt][j]>m:
                m = count[nt][j]
        rest = x - m
        scr.append(rest)
    return sum(scr)

####6 Pr, Prob calculator: $$needs profile

def Pr(Text, Profile):
    lst = []
    fin = 1
    n = len(Text)
    for i in range(n):
        lst.append(Profile[Text[i]][i])
        fin *= lst[i]
    return fin

####7 PMPP: $$needs profile and pr

def ProfileMostProbableKmer(Text,k,Profile):
    n = len(Text)
    mx = 0
    pos = 0
    for i in range(0,n-k+1):
        prob = Pr(Text[i:i+k],Profile)
        if prob > mx:
            mx = prob
            pos = i
    return Text[pos:pos+k]

####8 Greedy Motif Search: $$needs all above (except entropy)

def GreedyMotifSearch(Dna,k,t):
    n = len(Dna[0])
    best_motif = []
    for i in range(t):
        best_motif.append(Dna[i][0:k])
    for x in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][x:x+k])
        for y in range(1,t):
            p = Profile(Motifs[0:y])
            Motifs.append(ProfileMostProbableKmer(Dna[y], k, p))

        if Score(Motifs) < Score(best_motif):
            best_motif = Motifs
    return best_motif
