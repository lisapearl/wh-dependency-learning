# wh-dependency-learning
 For learning constraints on wh-dependencies. 
 Allows user to specify input files containing wh-dependency counts and 
 inital probabilities, along with the learning period, 
 how often learned probabilities are output, and
 what the smoothing constant is for relativizing probabilities.

If using, please cite the following:
	
Pearl, L. & Sprouse, J. 2013. 
Syntactic islands and learning biases: 
Combining experimental syntax and computational modeling 
to investigate the language acquisition problem. 
*Language Acquisition, 20, 23-68.* 
DOI 10.1080/10489223.2012.738742.


# Main code: online_learner.pl
  invoking learner with default options
  
  *online_learner.pl --inputfile <input-file-name> --chainfile <chain-file-name>*
 
##  run_online_learner
   Example executable file showing calls of online_learner.pl with different input files

## Input files expected
### inputfile (sample: input/child-all.counts) ### 
*expected input format*

*chain<\tab>frequency*

IP	560 

IP-VP	3755

IP-VP-CP-IP	7

IP-VP-CP-IP-VP	107

### chainfile (sample: input/child-all.chains) ###
*expected chainfile format*

*chain<\tab>initial_prob*

IP-VP-CP-IP-NP	0.5

IP-VP-CP-IP-VP-NP	0.5

## Output files generated ##
### Sample file: output/child-all.run1 ###
A single file containing a variety of information about the learned probabilities, and related information.

   
  

