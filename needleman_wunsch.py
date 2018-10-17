import numpy as np
import sys, os
import argparse
from subprocess import Popen, PIPE




class needleman_wunsch:
    '''
        Class containing methods required in implementation of Needleman-Wunsch algorithm
    '''
    
    def get_alignment_matrix(self,sequence1, sequence2):
        '''
                Function that fills up the 2D-matrix used for creation of final alignments.
                For a normal needleman_wunsch, the penalty values used are:
                
                mismatches: -3
                match     :  1
                gap       : -2
                
                @input:    sequence1             -- input reference sequence
                @input:    sequence2             -- input query sequence
                @input:    matrix                -- optional matrix, to be used for anchored version of the algorithm
                
                @return:   Matrix                -- constructed 2D matrix after applying the needleman_wunsch algorithm
        '''
        num_cols=len(sequence2)+1
        num_rows=len(sequence1)+1
    
        Matrix=np.zeros((num_rows,num_cols))
    
        for x in range(num_cols):
            #initialize first row. This happens when gaps are introduced at the beginning of sequence 1 without consuming any character of sequence 2
            Matrix[0][x]=(x*-2)
        
        for x in range(num_rows):
            #initialize first column. This happens when gaps are introduced at the beginning of sequence 2 without consuming any character of sequence 1
            Matrix[x][0]=(x*-2)

        for row in range(1,num_rows):
            for col in range(1,num_cols):
            
                penalty=1 if sequence1[row-1]==sequence2[col-1] else -3   #using hardcoded penalty scores
                match=Matrix[row-1][col-1]+penalty
                delete=Matrix[row-1][col]-2                           #gap is introduced in sequence 1
                insert=Matrix[row,col-1]-2                            #gap is introduced in sequence 2
            
                Matrix[row][col]=max(match,insert,delete)             #select the step which maximizes the overall score

        #print(Matrix)
        #print("Alignment_score: "+str(Matrix[num_rows-1][num_cols-1]))
        alignment_score=Matrix[num_rows-1][num_cols-1]                #final alignment score
        return (Matrix,alignment_score)


            

    def get_aligned_sequences(self,reference, query, alignment_matrix):
        '''
            Creates aligned sequences from the filled up alignment matrix.
        
            @input:    reference             -- input reference sequence
            @input:    query                 -- input query sequence
            @input:    alignment_matrix      -- alignment matrix constructed by get_alignment_matrix()
        
            @return:   result_ref            -- reference sequence after alignment
            @return:   result_que            -- query sequence after alignment
        
        '''
    
        result_ref=""
        result_que=""
    
        i=len(reference) #keeps track of indices of reference
        j=len(query)     #keeps track of indices of query
    
        while i>0 or j>0:
            #loop till we reach the start of reference or query
            penalty=1 if reference[i-1]==query[j-1] else -3
            d=-2 #penalty for a gap
                
            if i>0 and j>0 and alignment_matrix[i][j]==alignment_matrix[i-1][j-1]+penalty:     #if the current characters in both the sequences are consumed
                result_ref=reference[i-1]+result_ref
                result_que=query[j-1]+result_que
                i=i-1
                j=j-1
            elif i>0 and alignment_matrix[i][j]==alignment_matrix[i-1][j]+d:                   #if the current character in reference is consumed and a gap is inserted in the query
                result_ref=reference[i-1]+result_ref
                result_que="_"+result_que
                i=i-1
            else:                                                                               #if the current character in query sequence is consumed and a gap is inserted in the reference
                result_ref="_"+result_ref
                result_que=query[j-1]+result_que
                j=j-1

        return (result_ref,result_que)

    
    

    def process_input(self,filename):
        '''
                Parses input file and returns the input sequence
                @input:    reference             -- input reference sequence

                @return:   result_ref            -- reference sequence after alignment

                
        '''
        lines=[line.rstrip('\n') for line in open(filename)]
        #print(''.join(lines[1:]))
        return ''.join(lines[1:])
    


    def print_matrix(self,matrix):
        '''
                Prints filled up matrix
                @input:    matrix            -- input matrix
                
        '''
        for row in range(0,len(matrix)):
            print(matrix[row])




class driver:
    
    @classmethod
    def run_needleman_wunsch(cls,input_human,input_fly):
        '''
            Driver method to run needleman-wunsch alignment algorithm
        '''
        NMW=needleman_wunsch()
        reference_sequence=NMW.process_input(input_fly)
        #reference_sequence="ABCDE"
        #print("Reference sequence: \n"+reference_sequence+"\n")
    
        query_sequence=NMW.process_input(input_human)
        #query_sequence="ZABCDF"
        #print("Query sequence: \n"+query_sequence+"\n")
    
    
    
        M,alignment_score=NMW.get_alignment_matrix(reference_sequence,query_sequence)
    
        result_ref,result_que=NMW.get_aligned_sequences(reference_sequence,query_sequence,M)
    
    
    
        print('\n***************************************************\n')
        print("Running needleman wunsch: \n\n\n")
        print("Alignment score is:"+str(alignment_score))
        print("\n\nAligned reference is:\n"+result_ref)
        print("\n\nAligned Query is:\n"+result_que)
        print('\n***************************************************\n')






#driver.run_needleman_wunsch()


def make_arg_parser():
    parser = argparse.ArgumentParser(prog='needleman_wunsch.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-q","--query",
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to query fasta [required]")              #human sequence
    parser.add_argument("-r","--ref",
                        default=argparse.SUPPRESS,
                        required=True,
                        help="Path to reference fasta [required]")          #fly    sequence


    return parser




if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()


    driver.run_needleman_wunsch(args.query,args.ref)           # "HUMAN.PAX" "FLY.PAX"




