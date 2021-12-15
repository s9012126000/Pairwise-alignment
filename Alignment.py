import numpy as np
class sequence(object):
    '''
    Create object for storing DNA sequence
    Attribute
        seq: str for DNA sequence (eg. 'ATCG')
    '''
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
    '''
    Create object for storing score sets
    Attributes
        match: score for match DNA sequence
        mis: penalty for mismatch DNA sequence
        gap: penalty for DNA sequence Shifting
    '''
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
    '''
    Load 2 compared sequence object and build a score matrix for per-base consistency
    Attributes
        seq1: sequence object 1 queue to be align
        seq2: sequence object 2 queue to be align
        score: score object
    '''
    def __init__(self,seq1,seq2,score):
        self.seq1 = seq1
        self.seq2 = seq2
        self.score = score

    def score_Matrix(self):
        '''
        Create score matrix for per-base consistency by score object
        Return: a numpy array
        '''
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
    '''
    A global alignment aligns two sequences from beginning to end base on
    Needleman-Wunsch algorithm,aligning each letter in each sequence only
    once.An alignment is produced,regardless of whether or not there is
    similarity between the sequences.
    ref: NCBI
    '''
    def Global_start(self):
        '''
        a stepwise function call for global alignment
        '''
        score_arr = alignment.score_Matrix(self)
        main_arr = Global.Matrix(self, score_arr)
        self.align1, self.align2, self.link, self.value = Global.traceback(self, main_arr, score_arr)
        return Global.__str__(self)

    def Matrix(self,score_arr):
        '''
        Create the main matrix of largest accumulative score in 3 possible directions
        Args:
            score_arr: score matrix
        return:
            main matrix: an array
        '''
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
        '''
        return an array with default gap penalty
        '''
        for i in range(len(self.seq1) + 1):
            main_arr[i][0] = i * self.score.g()
        for j in range(len(self.seq2) + 1):
            main_arr[0][j] = j * self.score.g()
        return main_arr

    def traceback(self, main_arr, score_arr):
        '''
        Search and save the optimised path traceback from main matrix
        main_arr: main matrix
        Args:
            score_arr: score matrix
        return:
            alignment1: str, result of align1
            alignment2: str, result of align2
            link: str, a set of '|' indicate the location of match
            value: int, the score for alignment
        '''

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
    '''
    The Smith–Waterman algorithm performs local sequence alignment; that is,
     for determining most similar regions between two strings of nucleic acid sequences
    '''
    def Local_start(self):
        '''
        a stepwise function call for Local alignment
        '''
        score_arr = alignment.score_Matrix(self)
        main_arr, best, loc = Local.Matrix(self, score_arr)
        self.align1, self.align2, self.value, self.loc1, self.loc2 = Local.traceback(self, main_arr, score_arr, best, loc)
        return Local.__str__(self)

    def Matrix(self,score_arr):
        '''
        Create the main matrix of largest accumulative score in 3 possible directions,
        while the value <0, set the value into 0
        Args:
            score_arr: score matrix
        return:
            main matrix : an array
            best: int, the best value in the array
            loc: tuple, the location of best value in the array
        '''
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
        '''
        Search and save the optimised path of the most similar regions traceback from main matrix
        Args:
            main_arr: main matrix
            score_arr: score matrix
        return:
            alignment1: str, result of align1
            alignment2: str, result of align2
            value: int, the score for alignment
            aligned region1: tuple, coordinate for aligned region in alignment1
            aligned region2: tuple, coordinate for aligned region in alignment2
        '''
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
   
'''
Output results:

    Global alignment
        align 1 = ---T-CC-CAGT--TATGTCAGGGGACACGAGCATG-CAGAGAC
                     | || | ||  | | |||   | || |   ||| ||||  |
        align 2 = AATTGCCGCCGTCGT-TTTCA---G-CA-G-TTATGTCAGA-TC
        score = -36
        
    Local alignment
        Ref = CAGTTATGTCAG
        RefLoc: 4 -> 15
        seq = CAGTTATGTCAG
        SeqLoc: 22 -> 33
        score = 78
'''