import numpy as np

class sequence():
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

class Global(alignment):
    def Global_start(self):
        score_arr = alignment.score_Matrix(self)
        main_arr = Global.Matrix(self, score_arr)
        self.align1, self.align2, self.link, self.value = Global.traceback(self, main_arr, score_arr)
        return Global.__str__(self)

    def Matrix(self,score_arr):
        main_arr = np.zeros((len(self.seq1) + 1, len(self.seq2) + 1))
        main_arr01 = Global.gapPanalty(self, main_arr)
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                main_arr01[i][j] = max(main_arr01[i - 1][j - 1] + score_arr[i - 1][j - 1],
                                       main_arr01[i - 1][j] + self.score.g(),
                                       main_arr01[i][j - 1] + self.score.g()
                                       )
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

    def __str__(self):
        return '\n'+ \
               'align 1 = ' +str(self.align1)+'\n'+ \
               str(self.link)+'\n'+\
               'align 2 = '+ str(self.align2)+'\n'+\
               'score = ' + str(int(self.value))

class Local(alignment):
    def Local_start(self):
        score_arr = alignment.score_Matrix(self)
        main_arr, best, loc = Local.Matrix(self, score_arr)
        self.align1, self.align2, self.value, self.loc1, self.loc2 = Local.traceback(self, main_arr, score_arr, best, loc)
        return Local.__str__(self)

    def Matrix(self,score_arr):
        main_arr = np.zeros((len(self.seq1) + 1, len(self.seq2) + 1))
        best = 0
        loc = (0, 0)
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                main_arr[i][j] = max(0,main_arr[i - 1][j - 1] + score_arr[i - 1][j - 1],
                                       main_arr[i - 1][j] + self.score.g(),
                                       main_arr[i][j - 1] + self.score.g()
                                       )
                if main_arr[i][j] > best:
                    best = main_arr[i][j]
                    loc = (i, j)
        return main_arr, best, loc

    def traceback(self, main_arr, score_arr, best, loc):
        alignment1 = ''
        alignment2 = ''
        value = 0
        x = loc[0]
        y = loc[1]
        while best > 0:

            if main_arr[x][y] == main_arr[x - 1][y - 1] + score_arr[x - 1][y - 1]:
                alignment1 = self.seq1.sq()[x - 1] + alignment1
                alignment2 = self.seq2.sq()[y - 1] + alignment2
                value += main_arr[x][y]
                best = main_arr[x - 1][y - 1]
                x -= 1
                y -= 1

            elif main_arr[x][y] == main_arr[x - 1][y] + self.score.g():
                alignment1 = self.seq1.sq()[x - 1] + alignment1
                alignment2 = '-' + alignment2
                value += main_arr[x][y]
                best = main_arr[x - 1][y]
                x -= 1

            elif main_arr[x][y] == main_arr[x][y - 1] + self.score.g():
                alignment1 = '-' + alignment1
                alignment2 = self.seq2.sq()[y - 1] + alignment2
                value += main_arr[x][y]
                best = main_arr[x][y - 1]
                y -= 1

        return alignment1, alignment2, value, (x + 1, loc[0]), (y + 1, loc[1])

    def __str__(self):
        return '\n' + \
               'Ref = ' + str(self.align1) + '\n' + \
               'RefLoc: ' + str(self.loc1[0]) + ' -> ' + str(self.loc1[1]) + '\n' + \
               'seq = ' + str(self.align2) + '\n' + \
               'SeqLoc: '+ str(self.loc2[0])+' -> '+ str(self.loc2[1])+ '\n' + \
               'score = ' + str(int(self.value))

if __name__ == '__main__':
    seq1 = sequence("TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC")
    seq2 = sequence("AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC")
    score = score(1,-1,-1)
    Global_align = Global(seq1,seq2,score)
    Global_result = Global_align.Global_start()
    print(Global_result)
    Local_align = Local(seq1,seq2,score)
    Local_result = Local_align.Local_start()
    print(Local_result)