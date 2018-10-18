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





class anchored_needleman_wunsch(needleman_wunsch):
    '''
        Class containing methods required in implementation of anchored version of needleman wunsch algorithm
    '''

    def __init__(self):
        '''
            Constructor initializing the start and end indicies of matched regions in reference and query sequences
        '''
        self.indices_ref=[]
        self.indices_que=[]
    
    
    def set_indices(self,match_file):
        '''
            Extracts start and end indices of safe regions in both reference sequence (human protein sequence) and query sequence (fly protein sequence)
            
            @input: match_file     --source file containing information regarding the anchored regions
        '''
        lines=[line.rstrip('\n') for line in open(match_file)]
        
        for line in lines:
            values=line.split('\t')
            self.indices_ref.append([int(values[0]),int(values[1])])        #list of start,end indices of safe regions in human protein sequence
            self.indices_que.append([int(values[2]),int(values[3])])        #list of start,end indices of safe regions in fly protein sequence

    
    def get_aligned_sequences(self,reference, query, match_file):
        '''
            Performs an anchored version of the needleman wunsch algorithm
            @input:    reference             -- input reference sequence
            @input:    query                 -- input query sequence
            @input:    match_file            -- path to the match file, which contains anchor information
            
            @return:   result_ref            -- reference sequence after alignment
            @return:   result_que            -- query sequence after alignment
        
        '''
        self.set_indices(match_file)

        result_ref=""
        result_que=""
        curr_ref=0
        curr_que=0
        
        final_score=0
        
        for indices_ref,indices_que in zip(self.indices_ref,self.indices_que):
            #only consider the unsafe regions for alignment purposes
            start_ref=indices_ref[0]                                   #starting index for first anchor in input reference sequence
            end_ref=indices_ref[1]                                     #ending index for first anchor in input reference sequence
            start_que=indices_que[0]                                   #starting index for first anchor in input query sequence
            end_que=indices_que[1]                                     #ending index for first anchor in input reference sequence
            
            
            
            temp_ref=reference[curr_ref:start_ref]                      #substring of reference that requires alignment
            temp_que=query[curr_que:start_que]                          #substring of query that requires alignment

            M,score=super(anchored_needleman_wunsch,self).get_alignment_matrix(temp_ref,temp_que)               #compute alignment matrix
            
            final_score=final_score+score+(end_ref-start_ref+1)                                                         #assuming that the characters within the anchored regions are same in both sequences
            
            aligned_ref,aligned_que=super(anchored_needleman_wunsch,self).get_aligned_sequences(temp_ref,temp_que,M)    #compute alignment sequences for substrings

            result_ref=result_ref+aligned_ref+reference[start_ref:end_ref+1]                                            #add aligned substring and anchored substring to the result
            result_que=result_que+aligned_que+query[start_que:end_que+1]

            curr_ref=end_ref+1                                                                                          #set next curr_index for alignment
            curr_que=end_que+1


        if (curr_ref<=len(reference)-1) or (curr_que<=len(query)-1):        #align part of the reference and query sequences lying beyond the last anchor
            temp_ref=reference[curr_ref:]                      #substring of reference that requires alignment
            temp_que=query[curr_que:]                          #substring of query that requires alignment

            M,score=super(anchored_needleman_wunsch,self).get_alignment_matrix(temp_ref,temp_que)
            aligned_ref,aligned_que=super(anchored_needleman_wunsch,self).get_aligned_sequences(temp_ref,temp_que,M)
            
            final_score=final_score+score

            result_ref=result_ref+aligned_ref
            result_que=result_que+aligned_que
    


        return (result_ref,result_que,final_score)                                      #return final resulting aligned sequences along with final score






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



    @classmethod
    def run_needleman_wunsch_anchored(cls,input_human,input_fly,match_file):
        '''
                Driver method to run anchored version of needleman-wunsch alignment algorithm
        '''
        
        ANMW=anchored_needleman_wunsch()
        reference_sequence=ANMW.process_input(input_human)                   #to prevent out-of bounds error, human protein sequences are taken as reference in this algorithm
        #reference_sequence="ABCDE"
        #print("Reference sequence: \n"+reference_sequence+"\n")
        
        query_sequence=ANMW.process_input(input_fly)
        #query_sequence="ZABCDF"
        #print("Query sequence: \n"+query_sequence+"\n")

        
        result_ref,result_que, final_score=ANMW.get_aligned_sequences(reference_sequence,query_sequence,match_file)
        
        print('\n***************************************************\n')
        print("Running anchored needleman wunsch: \n\n\n")
        print("\n\nAligned reference is:\n"+result_ref)
        print("\n\nAligned Query is:\n"+result_que)
        print("\n\nFinal Score is:\n"+str(final_score))
        print('\n***************************************************\n')



    @classmethod
    def shuffle_and_run(cls):

        run_count=10000                                     #shuffle and run the implemented algorithm 10000 times.

        NMW=needleman_wunsch()
        reference_sequence=NMW.process_input("Fly_HOX.fa")
        query_sequence=NMW.process_input("Human_HOX.fa")
        M,original_alignment_score=NMW.get_alignment_matrix(reference_sequence,query_sequence)

        score_list=[]
        
        output_file=open('your_file.txt', 'w')


        import random

        while run_count:
            l_ref=list(reference_sequence)
            l_query=list(query_sequence)
            random.shuffle(l_ref)
            random.shuffle(l_query)
            shuffled_ref=''.join(l_ref)
            shuffled_query=''.join(l_query)

            shuffled_M,shuffled_score=NMW.get_alignment_matrix(shuffled_ref,shuffled_query)

            score_list.append(shuffled_score)
            
            print("run count: "+str(run_count)+" score: "+str(shuffled_score))

            run_count=run_count-1
            
            
            output_file.write("%s\n" % shuffled_score)

        output_file.close()
        



#driver.run_needleman_wunsch_anchored()
#driver.run_needleman_wunsch()

#driver.shuffle_and_run()

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
    parser.add_argument("-m","--match",
                        default=None,
                        required=False,
                        help="Path to match file [optional]")

    return parser




if __name__ == '__main__':
    parser = make_arg_parser()
    args = parser.parse_args()

    if args.match is None:
        driver.run_needleman_wunsch(args.query,args.ref)           # "HUMAN.PAX" "FLY.PAX"
    else:
        driver.run_needleman_wunsch_anchored(args.query,args.ref,args.match)



