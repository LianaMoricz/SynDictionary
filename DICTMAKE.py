from Bio.Seq import Seq

def make(amplicon, aftercutsite):
    dictlist = []
    l1 = len(amplicon)
    l2 = len(aftercutsite)
    l3 = l1 - l2

    lengthofeachitem = l3 // 2
    #if l3 % 2 != 0:
        #print("Delete base.")
        #return []
    
    #wholemutseq = amplicon[0:l3]

    seqprev = amplicon[0:lengthofeachitem] 
    seqafter = amplicon[lengthofeachitem:l3]


    for i in range(1, lengthofeachitem):

        newstring = seqafter[0:len(seqafter)-i]


        newstring = seqprev[-i:] + newstring

        newstringfin = newstring + aftercutsite
        dictlist.append(newstringfin)

    return dictlist


def main():
    # 5 to 3'
    amplicon = ("GCGAACTTTTGGCCGTGATGGGCAGTTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCGCAGGGCATCCAAGTATCGCCATCCGGGACTGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGG")
    #reference sequence around cut site (this stays at the end of all variants)
    aftercutsite = ("GCAGGGCATCCAAGTATCGCCATCCGGGAcTGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGG")

    #156 - 80 = 76
    #76 // 2 = 38 - 1 (non edited sequence) = 37

    dictlistfin = make(amplicon, aftercutsite)

    print("final list:")
    print()
    print("length =  " + str(len(dictlistfin)))
    for i, seq in enumerate(dictlistfin, 1):
    
        print(f"{seq}")


if __name__ == "__main__":
    main()








# amplicon = ("GCGAACTTTTGGCCGTGATGGGCAGTTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCGCAGGGCATCCAAGTATCGCCATCCGGGATGCGACTGCTCAATGGCCAACCTGTGGACGCCAAGGAGATGCAGGCCAGG")

# bparoundcs101 = ("TCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCTCCGGTGCCGGAAAGACGACCCTGCTGAATGCCCTTGCCTTTCGATCGCCG")
