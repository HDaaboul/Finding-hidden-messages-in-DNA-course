"""I put sorted programs from part I bioinformatics course (Stepik):
1-pattern counter: count occurances of a pattern in text
2-frequency mapping function: returns a dict with all k-mers with their frequency in text
3-most frequent: most appearing string pattern(s) of length k in text $note: needs freq map
4-reverse complement of DNA: return complement strand
5-pattern matching positioner: renders position(s) of input string in text
6-symbol array (two alg): counts occurances of a pattern in a window frame and records this with every step
(here we put the extended genome, genome + first half of genome and make the window as half genome)
the faster version counts appearances in a window and then with every slide if the first base
(which now is out of frame) was a target, then we subtract 1, and if the most recent base in slide
(the last added) is a target, then add 1.
7-CG skew array: measures skew of G-C, so for each C subtract 1 and for each G add 1
8-minimum skew: returns the position in which the skew has a minimum (nadir).
9-Hamming distance: the number of differences btw two strings
10-Approximate mismatch positioner: renders positions of pattern appearances in text with maximum mismatches of d
11-Approximate mismatch counter: counts freq of pattern of maximum mismatch of d
"""
####1 pattern counter

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

####2 freq map

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq:
            freq[Pattern] += 1
    return freq

####3 most frequent word $$ needs freq map

def FrequentWords(Text, k):
    list1 = []
    freq = FrequencyMap(Text,k)
    list2 = freq.values()
    mx = max(list2)
    for i in freq:
        if freq[i] == mx:
            list1.append(i)
    return list1

####4 reverse complement of DNA

def Reverse(Pattern):
    rev = Pattern[::-1]
    return rev
def Complement(Pattern):
    comp = []
    for i in range(len(Pattern)):
        x = Pattern[i]
        if x == "A":
            comp.append("T")
        elif x == "T":
            comp.append("A")
        elif x == "C":
            comp.append("G")
        elif x == "G":
            comp.append("C")
        else:
            comp.append("?")
    return ''.join(comp)
def ReverseComplement(Pattern):
    x = Complement(Pattern)
    return Reverse(x)

####5 Pattern matching positioner

def PatternMatching(Pattern, Genome):
    n = len(Genome)
    k = len(Pattern)
    list = []
    for i in range(n-k+1):
        x = Genome[i:i+k]
        if x == Pattern:
            list.append(i)
    return list

####6 symbol array $$ needs pattern counter

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = counter(symbol, ExtendedGenome[i:i+(n//2)])
    return array
#### 6 faster algorithm $$super imporant idea

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#### 7 CG skew array

def SkewArray(Genome):
    Skew = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    Skew[0] = 0
    for i in range(0,n):
        if ExtendedGenome[i] == "G":
            Skew[i+1] = Skew[i]+1
        if ExtendedGenome[i] == "C":
            Skew[i+1] = Skew[i]-1
        if ExtendedGenome[i] == "A" or ExtendedGenome[i] == "T":
            Skew[i+1] = Skew[i]
    return Skew.values()

####8 minimum skew $$needs CG skew array

def MinimumSkew(Genome):
    n = len(Genome)
    val = list(SkewArray(Genome).values())
    key = list(SkewArray(Genome).keys())
    mn = min(val)
    pos = []
    for i in range(n):
        if val[i] == mn:
            pos.append(key[i])
    return pos

####9 Hamming distance

def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

####10 Approximate mismatch positioner $$needs Hamming distance

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

####11 Approximate mismatch counter $$needs Hamming dist, pattern counter

def ApproximatePatternCount(Pattern, Text, d):
    n = len(Text)
    k = len(Pattern)
    count = 0
    for i in range(n-k+1):
        if HammingDistance(Text[i:i+k],Pattern) <= d:
            count += 1
    return count

