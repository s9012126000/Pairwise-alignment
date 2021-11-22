import numpy as np

class sequence(object):
    def __init__(self,seq):
        if type(seq) == str:
            self.seq = seq
        else:
            raise ValueError
    def sq(self):
        return self.seq
    def __str__(self):
        return '<Sequence object: '+str(self.seq)+'>'
    def __len__(self):
        return len(self.seq)

class score(object):
    def __init__(self,match,mis,gap):
        self.match = match
        self.mis = mis
        self.gap = gap
    def g(self):
        return self.gap
    def ms(self):
        return self.mis
    def ma(self):
        return self.match

class alignment(object):
    def __init__(self,seq1,seq2,score):
        self.seq1 = seq1
        self.seq2 = seq2
        self.score = score

    def startAlign(self):
        score_arr = alignment.score_Matrix(self)
        main_arr = alignment.Matrix(self,score_arr)
        align1, align2, link, value = alignment.traceback(self,main_arr, score_arr)
        return alignment.__str__(self, align1, align2, link, value)

    def score_Matrix(self):
        score_arr = np.zeros((len(self.seq1), len(self.seq2)))
        for i in range(len(self.seq1)):
            for j in range(len(self.seq2)):
                if self.seq1.sq()[i] == self.seq2.sq()[j]:
                    score_arr[i, j] = self.score.ma()
                else:
                    score_arr[i, j] = self.score.ms()
        self.score_arr = score_arr
        return self.score_arr

    def Matrix(self,score_arr):
        main_arr = np.zeros((len(self.seq1) + 1, len(self.seq2) + 1))
        main_arr01 = alignment.gapPanalty(self, main_arr)
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                main_arr01[i][j] = max(main_arr01[i - 1][j - 1] + score_arr[i - 1][j - 1],
                                       main_arr01[i - 1][j] + self.score.g(),
                                       main_arr01[i][j - 1] + self.score.g()
                                       )
        self.main_arr01 = main_arr01
        return main_arr01

    def gapPanalty(self,main_arr):
        for i in range(len(self.seq1) + 1):
            main_arr[i][0] = i * self.score.g()
        for j in range(len(self.seq2) + 1):
            main_arr[0][j] = j * self.score.g()
        return main_arr

    def traceback(self, main_arr, score_arr):
        alignment1 = ''
        alignment2 = ''
        link = ''
        refLenth = len(self.seq1)
        seqLenth = len(self.seq2)
        value = 0
        while refLenth > 0 or seqLenth > 0:

            if refLenth > 0 and seqLenth > 0 and main_arr[refLenth][seqLenth] == main_arr[refLenth-1][seqLenth-1] + score_arr[refLenth-1][seqLenth-1]:
                alignment1 = self.seq1.sq()[refLenth - 1] + alignment1
                alignment2 = self.seq2.sq()[seqLenth - 1] + alignment2
                if self.seq1.sq()[refLenth - 1]==self.seq2.sq()[seqLenth - 1]:
                    link = '|'+link
                else:
                    link = ' '+link
                value += main_arr[refLenth][seqLenth]
                refLenth -= 1
                seqLenth -= 1

            elif refLenth > 0 and main_arr[refLenth][seqLenth] == main_arr[refLenth - 1][seqLenth] + self.score.g():
                alignment1 = self.seq1.sq()[refLenth - 1] + alignment1
                alignment2 = '-' + alignment2
                link = ' ' + link
                value += main_arr[refLenth][seqLenth]
                refLenth -= 1

            else:
                alignment1 = '-' + alignment1
                alignment2 = self.seq2.sq()[seqLenth - 1] + alignment2
                link = ' ' + link
                value += main_arr[refLenth][seqLenth]
                seqLenth -= 1
        link = '          '+link
        return alignment1, alignment2, link, value

    def __str__(self,align1, align2, link, value):
        return '\n'+ \
               'align 1 = ' +align1+'\n'+ \
               link+'\n'+\
               'align 2 = '+ align2+'\n'+\
               'score = ' + str(int(value))


# if __name__ == '__main__':
#     seq1 = sequence("TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC")
#     seq2 = sequence("AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC")
#     score = score(1,-1,-1)
#     align = alignment(seq1,seq2,score)
#     result = align.startAlign()
#     print(result)
