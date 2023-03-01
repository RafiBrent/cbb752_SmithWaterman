
#!/usr/bin/python
__author__ = "Rafi Brent"
__email__ = "rafi.brent@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage: python3 hw1.py -i <input file> -s <score file>
### Example: python3 hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,
default=-2)
parser.add_argument('-e', '--extgap', help='extension gap',
required=False, default=-1)
args = parser.parse_args()

def runSW(inputFile, scoreFile, openGap=-2, extGap=-1):

  #read in the scoreFile to populate the simiilarity matrix and create a dictionary to associate residue ids with the correct index for the similarity matrix
  sim_mat=np.zeros((23,23), dtype=np.intc)
  index_dict={}
  with open(scoreFile) as f:
    lines = f.readlines()
    for i in range(len(lines)):
      if i !=0:
        temp_arr=lines[i].split()
        for j in range(len(temp_arr)):
          if j==0:
            index_dict[temp_arr[j]]=i-1
          else:
            sim_mat[i-1,j-1]=int(temp_arr[j])

  #read in the inputFile to obtain the two sequences to be aligned
  seq_1, seq_2='', ''
  with open(inputFile) as f:
    lines=f.readlines()
    seq_1=lines[0][:-1]
    seq_2=lines[1]
  if seq_2[-1]=='\n':
    seq_2=seq_2[:-1]

  #Sequentially fill a score matrix according to the Smith-Waterman algorithm, and save the path to "origin" location of each newly-created
  #score in a step matrix.
  score_mat=np.zeros((len(seq_1)+1,len(seq_2)+1), dtype=np.intc)
  steps_mat=np.zeros((len(seq_1)+1,len(seq_2)+1, 2), dtype=np.intc)
  for i in range(1, score_mat.shape[0]):
    for j in range(1, score_mat.shape[1]):

      #Zero is the minimum score possible in the Smith-Waterman algorithm
      max_val=0

      #Checking score obtained by aligning subsequent residues in each protein
      align_next=score_mat[i-1,j-1] + sim_mat[index_dict[seq_1[i-1]], index_dict[seq_2[j-1]]]
      if align_next>0:
        max_val=align_next
        steps_mat[i,j,0]=-1
        steps_mat[i,j,1]=-1

      #Checking scores obtained by inserting a gap in one of the sequences
      for k in range(1, i):
        if score_mat[i-k, j] + (k-1)*extGap +openGap >max_val:
          max_val=score_mat[i-k,j] +(k-1)*extGap +openGap
          steps_mat[i,j,0]=-k
          steps_mat[i,j,1]=0
      for c in range(1, j):
        if score_mat[i,j-c] +(c-1)*extGap +openGap >max_val:
          max_val=score_mat[i,j-c] + (c-1)*extGap +openGap
          steps_mat[i,j,0]=0
          steps_mat[i,j,1]=-c

      #Populate the current index of the score matrix with the highest possible score
      score_mat[i,j]=max_val


  #identify the location and value of the maximal alignment score
  max_ind=[0,0]
  max_score=0
  for i in range(1, score_mat.shape[0]):
    for j in range(1, score_mat.shape[1]):
      if score_mat[i,j]>=max_score:
        max_score=score_mat[i,j]
        max_ind=[i,j]


  #starting from the location of the maximal score, "retrace" the sequence alignment using the step matrix to insert appropriately-sized gaps where indicated.
  aligned_1, aligned_2 = '',''
  i=max_ind[0]
  j=max_ind[1]
  while True:
    step_i=steps_mat[i, j,0]
    step_j=steps_mat[i, j,1]

    #Retracing if the alignment originated without a gap
    if step_i==-1 and step_j==-1:
      aligned_1=seq_1[i-1] + aligned_1
      aligned_2=seq_2[j-1] + aligned_2
      j-=1
      i-=1

    #Retracing if the alignment originated with a gap
    elif step_i==0 and step_j<0:
      for ind in range(abs(step_j)):
        j-=1
        aligned_1= "-"+aligned_1
        aligned_2= seq_2[j]+aligned_2
    elif step_j==0 and step_i<0:
      for ind in range(abs(step_i)):
        i-=1
        aligned_2= "-"+aligned_2
        aligned_1= seq_1[i]+aligned_1
    
    #Terminate when the retracing is complete
    elif step_i==0 and step_j==0:
      break
    
    #This should never be triggered, but is included to avoid any possibility of an infinite loop
    else:
      print('something went wrong here')
      return

  #Add parentheses to indicate the "aligned" portions of the sequences, then add the "unaligned" portions of both sequences outside of these. 
  aligned_1 = '(' + aligned_1 + ')'
  aligned_2 = '(' + aligned_2 + ')'
  aligned_1=seq_1[0:i]+aligned_1
  aligned_2=seq_2[0:j]+aligned_2
  aligned_1=aligned_1 + seq_1[len(aligned_1.replace('-', ''))-2:]
  aligned_2=aligned_2 + seq_2[len(aligned_2.replace('-',''))-2:]

  #ensure that the "aligned" sequences are located at corresponding positions in both sequences after insertion of unaligned residues
  while aligned_1.find('(')>aligned_2.find('('):
    aligned_2 = ' ' + aligned_2
  while aligned_2.find('(')>aligned_1.find('('):
    aligned_1 = ' ' + aligned_1

  #include vertical bars to indicate identical aligned residues within the parentheses
  bars=''
  for i in range(min(len(aligned_1), len(aligned_2))):
    if aligned_1[i]==aligned_2[i] and i in range(aligned_1.find('(')+1, aligned_1.find(')')):
      bars+='|'
    else:
      bars+=' '

  #add spaces where necessary for compatibility with the sample output
  while len(aligned_1)<max(len(aligned_1), len(aligned_2)):
    aligned_1+=' '
  while len(bars)<max(len(aligned_1), len(aligned_2)):
    bars+=' '
  while len(aligned_2)<max(len(aligned_1), len(aligned_2)):
    aligned_2+=' '

  

  # Write the outputs to a text file named sw_output_rbrent.txt
  with open('sw_output_rbrent.txt', 'w') as f:

    #header including input sequences
    f.write('-----------')
    f.write('\n|Sequences|')
    f.write('\n-----------')
    f.write('\nsequence1')
    f.write('\n'+seq_1)
    f.write('\nsequence2')
    f.write('\n'+seq_2)
    f.write('\n--------------')
    f.write('\n|Score Matrix|')
    f.write('\n--------------')

    #Print a tab-delimited scoring matrix with associated row and column residue labels
    topline_str='\t\t'
    for res in seq_1:
      topline_str+=(res+'\t')
    f.write('\n'+topline_str + '\n')
    seq2_ind=0
    first_row=True
    for row in score_mat.transpose(): #transposed because I used seq_1 for rows and seq_2 for columns, the opposite of the sample outputs
      if first_row:
        first_row=False
        row_str='\t'
      else:
        row_str=seq_2[seq2_ind]+'\t'
        seq2_ind+=1
      for element in row:
        row_str+=(str(element)+'\t')
      f.write(row_str+'\n')

    #footer including final output
    f.write('----------------------')
    f.write('\n|Best Local Alignment|')
    f.write('\n----------------------')
    f.write('\nAlignment Score:' + str(max_score))
    f.write('\nAlignment Results:')
    f.write('\n'+aligned_1)
    f.write('\n'+bars+'\n')
    f.write(aligned_2+'\n')

#Running the function defined above with the user's arguments
runSW(args.input, args.score, float(args.opengap), float(args.extgap))
