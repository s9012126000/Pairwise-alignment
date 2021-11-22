import numpy as np


def score_Matrix(ref, seq, match, mismatch):
        score_arr = np.zeros((len(ref),len(seq)))
        for i in range(len(ref)):
            for j in range(len(seq)):
                if ref[i] == seq[j]:
                    score_arr[i, j] = match
                else:
                    score_arr[i, j] = mismatch
        return score_arr

def Matrix(score_arr,ref,seq,gap):
    main_arr = np.zeros((len(ref)+1, len(seq)+1))
    best = 0
    loc = (0,0)
    for i in range(1,len(ref)+1):
        for j in range(1,len(seq)+1):
            main_arr[i][j]= max(0,main_arr[i-1][j-1]+score_arr[i-1][j-1],
                                main_arr[i-1][j]+gap,
                                main_arr[i][j-1]+gap
                                )
            if main_arr[i][j]>best:
                best = main_arr[i][j]
                loc = (i,j)

    return main_arr, best, loc

def traceback(ref, seq, main_arr, score_arr, gap, score, loc):
    alignment1 = ''
    alignment2 = ''
    value = 0
    x = loc[0]
    y = loc[1]
    while score > 0:

        if main_arr[x][y] == main_arr[x-1][y-1]+score_arr[x-1][y-1]:
            alignment1 = ref[x - 1] + alignment1
            alignment2 = seq[y - 1] + alignment2
            value += main_arr[x][y]
            score = main_arr[x-1][y-1]
            x -= 1
            y -= 1
        elif main_arr[x][y] == main_arr[x-1][y] + gap:
            alignment1 = ref[x - 1] + alignment1
            alignment2 = '-'+ alignment2
            value += main_arr[x][y]
            score = main_arr[x - 1][y]
            x -= 1
        elif main_arr[x][y] == main_arr[x][y-1]+ gap:
            alignment1 = '-'+ alignment1
            alignment2 = seq[y - 1] + alignment2
            value += main_arr[x][y]
            score = main_arr[x][y-1]
            y -= 1

        if x > y:
            mergelen = loc[0]+len(ref)-loc[1]
            new_align1 = ref + (mergelen - len(ref)) * '-'
            new_align2 = (mergelen - len(seq))* '-' + seq
            link = ' '*(mergelen-x) + '|'*(loc[0]-x) + '' * (mergelen-loc[0])
        elif x == y:
            new_align1 = ref
            new_align2 = seq
        elif x < y:
            mergelen = loc[1]+len(seq)-loc[0]
            new_align1 = (mergelen - len(ref)) * '-' + ref
            new_align2 = seq + (mergelen - len(seq)) * '-'
            link = ' ' * (mergelen - loc[1]) + '|' * (loc[1] - y) + '' * (mergelen - y)
        elif len(ref)>len(seq):
            new_align1 = ref
            new_align2 = x * '-' + alignment2 + (len(ref)- loc[0]) * '-'
        elif len(ref)<len(seq):
            new_align1 = y * '-' + alignment1 + (len(seq)- loc[1]) * '-'
            new_align2 = seq

    return alignment1, alignment2, value, (x+1,loc[0]), new_align1, new_align2,link

if __name__ == '__main__':
    ref = "TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC"
    seq = "AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"
    match = 1
    mismatch = -1
    gap = -1
    score_arr = score_Matrix(ref,seq,match,mismatch)
    print(score_arr)
    main_arr, best, loc = Matrix(score_arr,ref,seq,gap)
    print(main_arr)
    align1,align2, value, location,a1, a2, link = traceback(ref, seq, main_arr, score_arr, gap, best, loc)
    print('')
    print('align 1 =',align1)
    print('align 2 =',align2)
    print('score = '+ str(int(value)))
    print(a1)
    print(link)
    print(a2)
    print(location)



    # seq1 = sequence("TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC")
    # seq2 = sequence("AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC")

