#!/usr/bin/perl

# online_learner.pl
# started on April 6, 2011 by Lisa Pearl

# goal: create an online Bayesian learning model that learns whether
# a particular movement sequence, as represented by a chain of 
# "container nodes", is grammatical.

# expected input format
# chain\tabfrequency
# Ex: from brown-adam.counts
# IP	560
# IP-VP	3755
# IP-VP-CP-IP	7
# IP-VP-CP-IP-VP	107
# IP-VP-CP-IP-VP-IP-VP	3

# expected chainfile format
# chain\tabinitial_prob
# Ex: from brown-adam.chains
# IP-VP-CP-IP-NP	0.5
# IP-VP-CP-IP-VP-NP	0.5


# invoking learner with default options
# online_learner.pl --inputfile <input-file-name> --chainfile <chain-file-name>

#################
# initialize variable values to be used
$chain_separator = "-";
$learner_runs = 1;
$data_points_to_run = 1000;
$print_incremental = "yes";
$print_interval = 100;
$smoothing_constant = 0.5; # based on Lidstone's Law
$alpha = 1.0; # used for Bayesian update, which we don't do anymore
$beta = 1.0; # used for Bayesian update, which we don't do anymore

# get user-defined command line options #
use Getopt::Long;
GetOptions("help|h" => \$help, # print brief help message
	             "chain_separator:s" => \$chain_separator, # separates container nodes on a chain (default = -)
	             "learner_runs:i" => \$learner_runs, # number of learners to run (default =1)
	             "data_points_to_run:i" => \$data_points_to_run, # total learning length (default = 1000 data points)
	             "print_incremental:s" => \$print_incremental, # print results every so many data points (default = yes)
	             "print_interval:i" => \$print_interval, # if printing incrementally, print every how many? (default = 100) [also enters these into gen_probs over time for each chain]
	             "inputfile=s" => \$inputfile, # input file containing chains and frequencies
	             "chainfile=s" => \$chainfile, # input file containing chains whose probability we want to track
	             "smoothing_constant:f" => \$smoothing_constant, # smoothing constant used (default = 0.5)
	              "alpha:f" => \$alpha, # alpha prior for update equations (default = 1.0)
	              "beta:f" => \$beta); # beta prior for update equations (default = 1.0)

if($help){
  die("Implements an online trigram-based learner for learning grammatical container node chain sequences representing movement.
        Options:
        help, h  = print brief help message listing options available
	chain_separator = optional string variable for what separates container nodes on a chain (default = -)
	learner_runs = optional integer variable for number of learners to run (default =1)
	data_points_to_run = optional integer variable for total learning length (default = 1000 data points)
	print_incremental:s = if set to yes, print results every so many data points (default = yes)
	print_interval = optional integer variable: if printing incrementally, print every how many? (default = 100) [also enters these into gen_probs over time for each chain]
	inputfile = mandatory string variable indicating input file containing chains and frequencies
	chainfile = mandatory string variable indicating input file containing chains whose probability we want to track
	smoothing_constant = optional float variable for smoothing constant used (default = 0.5)
	alpha = optional float variable for alpha prior for Bayesian update equations, if used (default = 1.0)
	beta = optional float variable for beta prior for Bayesian update equations, if used (default = 1.0)
");
}

# print out current variable values being used for this learning simulation
print("Options selected:
print_incremental = $print_incremental
print_interval = $print_interval
inputfile = $inputfile
chainfile = $chainfile
chain_separator = $chain_separator
learner_runs = $learner_runs
data_points_to_run = $data_points_to_run
smoothing_constant = $smoothing_constant
alpha = $alpha
beta = $beta\n");


#############################
##### main body of program #####
#############################
# (1) read in chain information from $inputfile
#       (a) initialize data structure that holds chains, based on chain nodes extracted
#            (should initialize with smoothing constant, rather than all 0s)
#       (b) print out list of chain nodes extracted from chains
# (2) read in chains to track from $chainfile
#       (a) initialize data structure that contains probabilities of chains to track
my ($ref_arr_learner_input, $ref_hash_node_types, $ref_hash_chains_to_track) = read_chain_info($inputfile, $chainfile, $chain_separator);
print_nodes($ref_hash_node_types);
my ($ref_trigrams, $total_trigrams, $total_tri_observations) = initialize_trigrams($ref_hash_node_types, $smoothing_constant);
my %observed_chains; # used in updating probabilities
my $num_obs_chains = 0; # used in updating probabilities


# create data structures that can more easily be passed to functions
$ref_print_options = {
		 print_incremental => $print_incremental,
		 print_interval => $print_interval
};

$ref_update_eq_vars = {
		  smoothing_constant => $smoothing_constant,
		  alpha => $alpha,
		  beta => $beta,
		  total_trigrams => $total_trigrams,
		  total_tri_observations => $total_tri_observations,
		  observed_chains => \%observed_chains,
		  num_obs_chains => $num_obs_chains
};


# (3) for data_points_to_run data points
#      (a) randomly grab a data point d, based on distribution from $inputfile
#      (b) update counts of chains, based on observation
#      (c) update probabilities of chains to track, based on d
#      (d) if printing_incremental, check to see if at print_interval
my $curr_dp;

for(my $run = 0; $run < $learner_runs; $run++){
  for(my $dp = 0; $dp < $data_points_to_run; $dp++){
    if($dp % 1000 == 0){ print(STDERR ".");}

    if(($print_incremental eq "yes") && ($dp % $print_interval == 0)){
      #print("$dp:\t");
      #print_chain_probs($ref_hash_chains_to_track);  # Bayesian learner probabilities
      print_chain_gen_probs($ref_hash_chains_to_track, $ref_trigrams, $ref_update_eq_vars, $dp); # trigram gen probabilities
    }

    $curr_dp = get_rand_dp($ref_arr_learner_input); # this is a string variable, like "IP-VP"
    #print("\n***debug: current data point pulled is $curr_dp\n");
    $ref_update_eq_vars->{"total_tri_observations"} = update_trigram_counts($ref_trigrams, $curr_dp, $chain_separator, $ref_update_eq_vars->{"total_tri_observations"});
    #print("debug: total trigram observations at this point = $ref_update_eq_vars->{\"total_tri_observations\"}\n");
    update_observed_chains($ref_update_eq_vars, $curr_dp);
    # ****if trying to do a Bayesian learner***
    #update_chain_probabilities($ref_hash_chains_to_track, $ref_trigams, $ref_update_eq_vars, $curr_dp);
    #print("\n***debug: next data point\n");
  }
  # (4) print out final probability results
  # for each chain in chains_to_track, print out p_grammatical
  #print("$data_points_to_run:\t");
  #print_chain_probs($ref_hash_chains_to_track);
  #print("debug: full chain info\n");
  #print_full_chain_info($ref_hash_chains_to_track);
  print("debug: chain gen probabilities\n");
  print_chain_gen_probs($ref_hash_chains_to_track, $ref_trigrams, $ref_update_eq_vars, $data_points_to_run);
  print("debug: chain gen probabilities over time\n");
  print_all_chain_gen_probs($ref_hash_chains_to_track);
  print("debug: observed data\n");
  print_hashref($ref_update_eq_vars->{"observed_chains"});
  print("debug: final trigram counts\n");
  print_3Dhashref($ref_trigrams);
  
}

##########################
###### learner functions #####
##########################

sub update_chain_probabilities{
  my ($ref_hash_chains_to_track, $ref_trigams, $ref_update_eq_vars, $curr_dp) = @_;
  my $ref_chain_props;
  # for each chain to track, we update p_grammatical, based on trigrams observed so far
  foreach my $chain (sort keys (%$ref_hash_chains_to_track)){
    #print("\n***\ndebug: updating probability for chain $chain\n");
    update_chain_prob($ref_hash_chains_to_track->{$chain}, $ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain);
    #print("debug: p_grammatical for $chain = $ref_hash_chains_to_track->{$chain}->{\"p_grammatical\"}\n");
  }
}

sub update_chain_prob{
  my ($ref_hash_chain, $ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain) = @_;

  # general form of updating p_grammatical
  # p_grammatical = (alpha + data_x)/(alpha + beta + totaldata_x)
  #
  # totaldata_x is always updated by 1 for having seen one chain data point
  $ref_hash_chain->{"totaldata_x"} += 1;
  # data_x is updated by phi_x = prob(chain = grammatical | observed data) = p(g | o)
  #   we calculate phi_x using Bayes rule (likelihood * prior / evidence)
  #   = [p(o | g) * p(g)] / [p(o|g)*p(g) + p(o|not g) * p(not g)]
  my $phi_x = calculate_phi_x($ref_hash_chain->{"p_grammatical"}, $ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain);
  #print("debug: phi_x is $phi_x\n");

  $ref_hash_chain->{"data_x"} += $phi_x;

  # now calculate new p_grammatical
  my $a = $ref_update_eq_vars->{"alpha"};
  my $b = $ref_update_eq_vars->{"beta"};
  my $d_x = $ref_hash_chain->{"data_x"};
  my $td_x = $ref_hash_chain->{"totaldata_x"};
  #print("debug, inside update_chain_prob: calculating ($a + $d_x)/($a + $b + $td_x)\n");
  $ref_hash_chain->{"p_grammatical"} = ($a + $d_x)/($a + $b + $td_x);

  #print("debug: properties of current chain after update\n");
  #print("debug: data_x = $ref_hash_chain->{\"data_x\"}\n");
  #print("debug: totaldata_x = $ref_hash_chain->{\"totaldata_x\"}\n");
  #print("debug: p_grammatical = $ref_hash_chain->{\"p_grammatical\"}\n");
}

sub calculate_phi_x{
  my($p_grammatical, $ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain) = @_;

  # calculating phi_x = [p(o | g) * p(g)] / [p(o|g)*p(g) + p(o|not g) * p(not g)]
  my $phi_x;

  # prior p(g) = $p_grammatical
  my $prior = $p_grammatical;
  #print("debug: p_grammatical = $p_grammatical\n");

  # likelihood p(o|g) = likelihood of observed data given that chain is grammatical
  # because probabilities of generating chain are so small, doing a log transform
  #   more specifically, -1/log(generating the chain from trigrams)
  #   applied to a weighted hypothesis space of observed chains
  #   i.e., calculating the likelihood of observing the current chain
  #          given the set of chains that could have been observed, 
  #          given how easy/hard to generate them
  #print("debug: calculating likelihood of observed data point $curr_dp given that $chain is grammatical\n");
  my $likelihood = calculate_likelihood($ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain, 1);
  #print("debug: numerator likelihood = $likelihood\n");

  # evidence = p(o|g)*p(g) + p(o|not g) * p(not g)
  #                = likelihood*prior [chain is grammatical] + [what happens when chain is not grammatical]
  my $evidence = $prior * $likelihood; # first term where chain is grammatical
  #print("debug: evidence begins initialized to numerator (prior * likelihood) = $prior * $likelihood = $evidence\n");
  # if curr_dp = chain in question, p(o | not g) = 0 since chain has been generated so it must be grammatical
  # so evidence is not updated any further
  # otherwise, add in  (1-prior)* adjusted likelihood calculation that doesn't include chain in question
  if($curr_dp ne $chain){
    #print("debug: $curr_dp is not $chain, so calculating likelihood if $chain not allowed\n");
    $likelihood_ungram = calculate_likelihood($ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain, 0);
    #print("debug: $curr_dp is not $chain, so adding (1-$prior)*($likelihood_ungram) to $evidence\n");
    $evidence += (1-$prior)*($likelihood_ungram);
  }
  #print("debug: evidence = $evidence\n");

  $phi_x = ($prior*$likelihood)/$evidence;
  #print("debug: calculated phi_x = ($prior*$likelihood)/$evidence = $phi_x\n");
  
  return $phi_x;
}

sub calculate_likelihood{
  my($ref_trigrams, $ref_update_eq_vars, $curr_dp, $chain, $with_chain) = @_;
  my $likelihood;

  # likelihood formula = prob(gen this data point)/total_prob(gen all observed data points + chain-in-question if not seen yet)
  # similar to Lidstone's Law?
  # [gen * (freq(curr_dp))]/ ([for all observed]gen * (freq(those_observed_dp)))  possibly + [chain_in_question](gen * smoothing_constant)

  # first part: probability of generating observed data point
  my $gen_curr_dp = calculate_gen_prob($ref_trigrams, $ref_update_eq_vars, $curr_dp);
  #print("debug: log_gen of current dp $curr_dp = $log_gen_curr_dp\n");
  #print("debug: gen_prob of current dp $curr_dp = $gen_curr_dp\n");

  # second part: observed frequency of curr_dp in observations
  my $freq_curr_dp = $ref_update_eq_vars->{"observed_chains"}->{$curr_dp};
  #print("debug: observed frequency of $curr_dp is $freq_curr_dp\n");

  my $numerator = $gen_curr_dp * $freq_curr_dp;
  #my $numerator = $freq_curr_dp;
  #print("debug: likelihood numerator for curr dp $curr_dp = $gen_curr_dp * $freq_curr_dp = $numerator\n");
  #print("debug: likelihood numerator for curr dp $curr_dp = $freq_curr_dp = $numerator\n");
  my $denominator = $numerator; # initialize with this value since this is one of the terms to be added in

  # third part: generating all observed chains
  # fourth part: frequency of all observed chains
  foreach my $obs_chain (sort (keys %{$ref_update_eq_vars->{"observed_chains"}})){
    #print("debug: observed chain being considered = $obs_chain\n");
    if($with_chain || 
       (!$with_chain && ($obs_chain ne $chain))){ # in case we're excluding the chain in question
      if($obs_chain ne $curr_dp){ # already added this term in (it's the same as the numerator)
	$gen_curr_dp = calculate_gen_prob($ref_trigrams, $ref_update_eq_vars, $obs_chain);
	#print("debug: gen_prob of $obs_chain = $gen_curr_dp\n");
	$freq_curr_dp = $ref_update_eq_vars->{"observed_chains"}->{$obs_chain};
	#print("debug: observed frequency of $obs_chain is $freq_curr_dp\n");
	$denominator +=  $gen_curr_dp * $freq_curr_dp;
	#$denominator += $freq_curr_dp;
	#print("debug: adding $gen_curr_dp * $freq_curr_dp for $obs_chain, new denom = $denominator\n");
	#print("debug: adding $freq_curr_dp for $obs_chain, new denom = $denominator\n");
      }
    }
  }

  #print("debug: calculating with chain in question $chain? $with_chain\n");
  # fifth part: optional addition of unobserved chain-in-question, since could still occur is calculating with_chain
  if($with_chain && !exists($ref_update_eq_vars->{"observed_chains"}->{$chain})){
    #print("debug: chain in question $chain has not been observed yet - adding default to denominator, curr val = $denominator\n");
    my $gen_chain = calculate_gen_prob($ref_trigrams, $ref_update_eq_vars, $chain);
    my $freq_chain = $ref_update_eq_vars->{"smoothing_constant"};
    $denominator +=  $gen_chain * $freq_chain;
    #$denominator += $freq_chain;
    #print("debug: adding $gen_chain * $freq_chain for $chain, new denom = $denominator\n");
    #print("debug: adding $freq_chain for $chain, new denom = $denominator\n");
  }

  $likelihood = $numerator/$denominator;
  #print("debug: likelihood = $numerator/$denominator = $likelihood\n");

  return $likelihood;
}

sub calculate_gen_prob{
  my ($ref_trigrams, $ref_update_eq_vars, $chain) = @_;
  my $gen_prob = 1;

  # calculate probability of generating chain, given sequence of trigrams
  # break down $chain into sequences of trigrams 
  #    and add their log(probabilities)
  my @chain_nodes = split(/$chain_separator/, $chain);
  
  # check to see if only one node long
  if($#chain_nodes == 0){
    # if so, just doing start-node-end
    # divide frequency by total trigram counts
    #print("debug, in calculate gen_prob for $chain: multiplying ($ref_trigrams->{\"start\"}->{$chain_nodes[0]}->{\"end\"}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $gen_prob *= ($ref_trigrams->{"start"}->{$chain_nodes[0]}->{"end"}/$ref_update_eq_vars->{"total_tri_observations"});
  }else{
    # do start-node0-node1 through node(i-1)-node(i)-end
    # do start
    #print("debug, calculating gen_prob for $chain: multiplying gen prob for start-$chain_nodes[0]-$chain_nodes[1]: ($ref_trigrams->{\"start\"}->{$chain_nodes[0]}->{$chain_nodes[1]}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $gen_prob *= ($ref_trigrams->{"start"}->{$chain_nodes[0]}->{$chain_nodes[1]}/$ref_update_eq_vars->{"total_tri_observations"});

    #all but end
    for(my $i = 0; $i < $#chain_nodes -1; $i++){
      #print("debug, calculating gen_prob for $chain: multiplying gen prob for $chain_nodes[$i]-$chain_nodes[$i+1]-$chain_nodes[$i+2]: ($ref_trigrams->{$chain_nodes[$i]}->{$chain_nodes[$i+1]}->{$chain_nodes[$i+2]}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");      
      $gen_prob *= ($ref_trigrams->{$chain_nodes[$i]}->{$chain_nodes[$i+1]}->{$chain_nodes[$i+2]}/$ref_update_eq_vars->{"total_tri_observations"});      
    }
    # do end
    #print("debug, calculating gen_prob for $chain: multiplying gen prob for $chain_nodes[$#chain_nodes-1]-$chain_nodes[$#chain_nodes]-end: ($ref_trigrams->{$chain_nodes[$#chain_nodes-1]}->{$chain_nodes[$#chain_nodes]}->{\"end\"}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $gen_prob *= ($ref_trigrams->{$chain_nodes[$#chain_nodes-1]}->{$chain_nodes[$#chain_nodes]}->{"end"}/$ref_update_eq_vars->{"total_tri_observations"});
  }

  #print("debug: gen_prob = $gen_prob\n");
  return $gen_prob;
}

# don't use this anymore --> basically does the log of the gen prob calculation
sub calculate_log_gen{
  my ($ref_trigrams, $ref_update_eq_vars, $chain) = @_;
  my $log_gen = 0;

  # calculate probability of generating chain, given sequence of trigrams
  # break down $chain into sequences of trigrams 
  #    and add their log(probabilities)
  my @chain_nodes = split(/$chain_separator/, $chain);
  
  # check to see if only one node long
  if($#chain_nodes == 0){
    # if so, just doing start-node-end
    # divide frequency by total trigram counts
    #print("debug, in calculate log_gen for $chain: calculating log($ref_trigrams->{\"start\"}->{$chain_nodes[0]}->{\"end\"}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $log_gen += log($ref_trigrams->{"start"}->{$chain_nodes[0]}->{"end"}/$ref_update_eq_vars->{"total_tri_observations"});
  }else{
    # do start-node0-node1 through node(i-1)-node(i)-end
    # do start
    #print("debug, calculating log_gen for $chain: adding gen prob for start-$chain_nodes[0]-$chain_nodes[1]: log($ref_trigrams->{\"start\"}->{$chain_nodes[0]}->{$chain_nodes[1]}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $log_gen += log($ref_trigrams->{"start"}->{$chain_nodes[0]}->{$chain_nodes[1]}/$ref_update_eq_vars->{"total_tri_observations"});

    #all but end
    for(my $i = 0; $i < $#chain_nodes -1; $i++){
      #print("debug, calculating log_gen for $chain: adding gen prob for $chain_nodes[$i]-$chain_nodes[$i+1]-$chain_nodes[$i+2]: log($ref_trigrams->{$chain_nodes[$i]}->{$chain_nodes[$i+1]}->{$chain_nodes[$i+2]}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");      
      $log_gen += log($ref_trigrams->{$chain_nodes[$i]}->{$chain_nodes[$i+1]}->{$chain_nodes[$i+2]}/$ref_update_eq_vars->{"total_tri_observations"});      
    }
    # do end
    #print("debug, calculating log_gen for $chain: adding gen prob for $chain_nodes[$#chain_nodes-1]-$chain_nodes[$#chain_nodes]-end: log($ref_trigrams->{$chain_nodes[$#chain_nodes-1]}->{$chain_nodes[$#chain_nodes]}->{\"end\"}/$ref_update_eq_vars->{\"total_tri_observations\"})\n");
    $log_gen += log($ref_trigrams->{$chain_nodes[$#chain_nodes-1]}->{$chain_nodes[$#chain_nodes]}->{"end"}/$ref_update_eq_vars->{"total_tri_observations"});
  }

  #print("debug: log_gen before transform for $chain = $log_gen\n");
  #$log_gen = 1/(-1*$log_gen);
  #print("debug: log_gen after transform for $chain = $log_gen\n");

  return $log_gen;
}

sub update_observed_chains{
  my ($ref_update_eq_vars, $curr_dp) = @_;
  # increment  %observed_chains and $num_obs_chains
  if(!exists($ref_update_eq_vars->{"observed_chains"}->{$curr_dp})){
    # smoothing constant added in for all observed chains
    $ref_update_eq_vars->{"observed_chains"}->{$curr_dp} = $ref_update_eq_vars->{"smoothing_constant"};
    $ref_update_eq_vars->{"num_obs_chains"} += $ref_update_eq_vars->{"smoothing_constant"};
  }
  $ref_update_eq_vars->{"observed_chains"}->{$curr_dp} += 1;
  $ref_update_eq_vars->{"num_obs_chains"} += 1;
  #print("debug, after updating observed chains:\n");
  #print("debug: observed_chains = \n");
  #print_hashref($ref_update_eq_vars->{"observed_chains"});
  #print("debug: num_obs_chains = $ref_update_eq_vars->{\"num_obs_chains\"}\n");
}

sub update_trigram_counts{
  my ($ref_trigrams, $curr_dp, $chain_separator, $total_obs) = @_;

  # break down $curr_dp into sequences of trigrams & add 1 to all counts
  my @chain_nodes = split(/$chain_separator/, $curr_dp);
  # for each sequence of 3, update trigrams
  for(my $i = 0; $i < $#chain_nodes-1; $i++){
    #print("debug, trigram update: $chain_nodes[$i], $chain_nodes[$i+1], $chain_nodes[$i+2]\n");
    $ref_trigrams->{$chain_nodes[$i]}->{$chain_nodes[$i+1]}->{$chain_nodes[$i+2]} += 1;
    $total_obs += 1;
  }
  # also include sequence start-node1-node2 (where node2 might be end)
  # and include sequence node(i-1)-node(i)-end (where node(i-1) might be start)
  # check to see if chain is 1 or 2 nodes long
  if($#chain_nodes == 0){
    # just one node so update start-node0-end
    $ref_trigrams->{"start"}->{$chain_nodes[0]}->{"end"} += 1;
    #print("debug, trigram update: start, $chain_nodes[0], end\n");
    $total_obs += 1;
  }elsif($#chain_nodes == 1){
    # two nodes so update start-node0-node1 and node0-node1-end
    $ref_trigrams->{"start"}->{$chain_nodes[0]}->{$chain_nodes[1]} += 1; 
    #print("debug, trigram update: start, $chain_nodes[0], $chain_nodes[1]\n");
    $ref_trigrams->{$chain_nodes[0]}->{$chain_nodes[1]}->{"end"} += 1; 
    #print("debug, trigram update: $chain_nodes[0], $chain_nodes[1], end\n");
    $total_obs += 2;
  }else{
    # normal update with start-node0-node1 and node(i-1)-node(i)-end
    $ref_trigrams->{"start"}->{$chain_nodes[0]}->{$chain_nodes[1]} += 1; 
    #print("debug, trigram update: start, $chain_nodes[0], $chain_nodes[1]\n");
    $ref_trigrams->{$chain_nodes[$#chain_nodes-1]}->{$chain_nodes[$#chain_nodes]}->{"end"} += 1; 
    #print("debug, trigram update: start, $chain_nodes[$#chain_nodes-1], $chain_nodes[$#chain_nodes], end\n");
    $total_obs += 2;
  }
  #print("debug: current state of trigrams\n");
  #print_3Dhashref($ref_trigrams);
  return $total_obs;
}

sub get_rand_dp{
  my($ref_arr_learner_input) = @_;
  my $curr_dp;

  # should pull out an index of a data point
  my $dp_index = int(rand($#{$ref_arr_learner_input}+1));
  #print("debug: index for random data point is $dp_index\n");

  return $ref_arr_learner_input->[$dp_index];
}


##########################
### initialization functions ###
#########################

sub initialize_trigrams{
  my ($ref_hash_node_types, $smoothing_constant) = @_;

  # trigram data structure will be a 3D hash {node1}{node2}{node3} = frequency
  # all the entries will be initialized to the smoothing_constant

  # get list of node types (don't care how many there were, really)
  my @node_types = keys(%$ref_hash_node_types);
  my $num_node_types = $#node_types + 1; # last index + 1

  # get total number of trigrams (will be useful for learning calculations later on)
  my $total_trigrams = 0;
  my $total_obs = 0; # running total of observations

  # initialize, remembering to include trigrams for start and end
  my %trigrams;

  # do main body
  foreach my $node1 (@node_types){
    foreach my $node2 (@node_types){
      foreach my $node3 (@node_types){
	$trigrams{$node1}{$node2}{$node3} = $smoothing_constant;
	$total_obs += $smoothing_constant;
	$total_trigrams += 1;
      }
    }
  }

  # do start
      foreach my $node2 (@node_types){
      foreach my $node3 (@node_types){
	$trigrams{"start"}{$node2}{$node3} = $smoothing_constant;
	$total_obs += $smoothing_constant;
	$total_trigrams += 1;
      }
    }

  # do end
      foreach my $node1 (@node_types){
      foreach my $node2 (@node_types){
	$trigrams{$node1}{$node2}{"end"} = $smoothing_constant;
	$total_obs += $smoothing_constant;
	$total_trigrams += 1;
      }
    }

  # do singletons
  foreach my $node1 (@node_types){
	$trigrams{"start"}{$node1}{"end"} = $smoothing_constant;
	$total_obs += $smoothing_constant;
	$total_trigrams += 1;
      }


  #debug: make sure trigrams hash contains what we think it does
  #print_3Dhashref(\%trigrams);

  return(\%trigrams, $total_trigrams, $total_obs);
}

######################
### input functions ######
######################

sub read_chain_info{
  my ($inputfile, $chainfile, $chain_separator) = @_;

  # first read in chains and frequencies
  # put into learner_input data structure (array, where index i = chain, based on frequency)
  my ($ref_arr_learner_input, $ref_hash_node_types) = read_chain_freq($inputfile, $chain_separator);

  # then read in chains to track, and their initial probabilities of being grammatical
  my %chains_to_track = read_chains_to_track($chainfile, $chain_separator, $ref_hash_node_types);

  return ($ref_arr_learner_input, $ref_hash_node_types, \%chains_to_track);
}

sub read_chains_to_track{
  my ($chainfile, $chain_separator, $ref_hash_node_types) = @_;

  # this involves initializing a hash of hash references that contain the various quantities needed in the update equations:
  # (1) data_x [float: number of observations of this chain]
  # (2) totaldata_x [float: number of informative observations with respect to this chain]
  # (3) p_grammatical [float: current probability that chain is grammatical]
  # (4) gen probs: hash that tracks probability of generating chain, based on number of data points seen
  my %chains_to_track = ();
  my $chain, $init_prob;
  open(CHAINS, "$chainfile") || die("Couldn't open $chainfile that contains chains to track and their initial probabilities of being grammatical.\n");
  while(defined($chainline = <CHAINS>)){
    if($chainline =~ /(\S+)\s+(\d+\.?\d*)/){
      $chain = $1;
      # extract any additional chain nodes that weren't seen in the input
      extract_chain_nodes($ref_hash_node_types, $chain, $chain_separator);
      $init_prob = $2;
      # create anonymous hashref entry containing relevant entries
      my %gen_probs = ();

      $chains_to_track{$chain} = {
				  data_x => 0,
				  totaldata_x => 0,
				  p_grammatical => $init_prob,
				  gen_probs => \%gen_probs,
				  };
    }
  }
  close(CHAINS);
  #debug
  #print_hash_of_hashrefs(%chains_to_track);
  return %chains_to_track;
}

sub read_chain_freq{
  my ($inputfile, $chain_separator) = @_;

  my @learner_input = ();
  my $curr_index = 0;
  my $chain, $freq;
  my %node_types = (); # stores container node types we have extracted from input sequences

  open(IN, "$inputfile") || die("Couldn't open $inputfile that contains input chains and frequencies.\n");
  while(defined($inline = <IN>)){
    if($inline =~ /(\S+)\s+(\d+)/){
      $chain = $1;
      $freq = $2;
      # extract container node types from chain in case there are new ones
      extract_chain_nodes(\%node_types, $chain, $chain_separator);

      # initialize @learner_input with chain frequency entries containing chain
      for($i = $curr_index; $i < $curr_index+$freq; $i++){
	$learner_input[$i] = $chain;
      }
      # update current index
      $curr_index = $curr_index+$freq;
    }
  }

  # debug check
  #print("debug: after reading in $inputfile, learner_input looks like this:\n");
  #for($i = 0; $i < $#learner_input; $i++){
  #  print("learner_input[$i] =  $learner_input[$i]\n");
  #}

  close(IN);
  return (\@learner_input, \%node_types);
}

sub extract_chain_nodes{
  my ($ref_hash_node_types, $chain, $chain_separator) = @_;
  my $node;
  # split current chain, and put into node_types hash
  my @chain_nodes = split(/$chain_separator/, $chain);
  foreach $node (@chain_nodes){
    # create if doesn't exist
    if(!exists($ref_hash_node_types->{$node})){
      #print("debug: new node type not encountered before = $node\n");
      $ref_hash_node_types->{$node} = 0;
    }
    # increment by 1
    $ref_hash_node_types->{$node}++
  }
}

######################
#### print functions #####
######################
sub print_full_chain_info{
  my ($ref_hash_chains_to_track) = @_;
    foreach my $chain (sort keys (%$ref_hash_chains_to_track)){
      print("$chain, p_grammatical: $ref_hash_chains_to_track->{$chain}->{\"p_grammatical\"}\t");
      print("$chain, data_x: $ref_hash_chains_to_track->{$chain}->{\"data_x\"}\t");      
      print("$chain, totaldata_x: $ref_hash_chains_to_track->{$chain}->{\"totaldata_x\"}\t");    
      print("\n");
    }
}

sub print_chain_gen_probs{
  my($ref_hash_chains_to_track, $ref_trigrams, $ref_update_eq_vars, $dp) = @_;
  my $gen_chain;
  print("$dp:\n");
  foreach my $chain (sort keys %$ref_hash_chains_to_track){
    $gen_chain = calculate_gen_prob($ref_trigrams, $ref_update_eq_vars, $chain);
    print("$chain gen probability = $gen_chain\n");
    # also update gen_probs for each chain
    $ref_hash_chains_to_track->{$chain}->{"gen_probs"}->{$dp} = $gen_chain;
    #print("debug: ref_hash_chains_to_track->{$chain}->{\"gen_probs\"}->{$dp} is now $ref_hash_chains_to_track->{$chain}->{\"gen_probs\"}->{$dp}\n");
  }
}

sub  print_all_chain_gen_probs{
  my ($ref_hash_chains_to_track) = @_;

  # print out initial listing of increments
  print_gen_prob_keys($ref_hash_chains_to_track);

  foreach my $chain (sort keys %$ref_hash_chains_to_track){
    print("$chain\t");
    #print("debug: now printing gen prob hashref for ref_hash_chains_to_track->{$chain}: $ref_hash_chains_to_track->{$chain}\n");
    print_hashref_num($ref_hash_chains_to_track->{$chain}->{"gen_probs"});
  }
}

sub print_gen_prob_keys{
  my ($ref_hash_chains_to_track) = @_;
  print("chain:\t");
  my @chains = keys %$ref_hash_chains_to_track;
  # just do for first chain available
  foreach my $key (sort {$a <=> $b} keys (%{$ref_hash_chains_to_track->{$chains[0]}->{"gen_probs"}})){
    print("$key\t");
  }
  print("\n");
}

sub print_hashref_num{
  my ($my_hashref) = @_;

  # note: :$a <=> $b yields a numerical sort by ascending order
  foreach $key (sort {$a <=> $b} keys (%$my_hashref)){
    print("$my_hashref->{$key}\t");
  }
  print("\n");
}

sub print_chain_probs{
  my ($ref_hash_chains_to_track) = @_;
  print("Current chain probabilities\n");
    foreach my $chain (sort keys (%$ref_hash_chains_to_track)){
      print("$chain: $ref_hash_chains_to_track->{$chain}->{\"p_grammatical\"}\t");
    }
  print("\n");
}

sub print_hash_of_hashrefs{
  my (%my_hash) = @_;

  while(($key, $hashref) = each(%my_hash)){
    print("$key:\n ");
    while(($refkey, $refval) = each(%$hashref)){
      print("$refkey = $refval\n");
    }
  }
}

sub print_hashref{
  my ($my_hashref) = @_;

  foreach $key (sort keys (%$my_hashref)){
    print("$key = $my_hashref->{$key}\n");
  }
}

sub print_3Dhashref{
  my ($ref_to_3Dhash) = @_;

  foreach my $key1 (sort keys (%$ref_to_3Dhash)){
    foreach my $key2 (sort keys (%{$ref_to_3Dhash->{$key1}})){
      foreach my $key3 (sort keys (%{$ref_to_3Dhash->{$key1}->{$key2}})){
	print("$key1 $key2 $key3 $ref_to_3Dhash->{$key1}->{$key2}->{$key3}\n");
      }
    }
  }
}

sub print_nodes{
  my ($ref_hash_node_types) = @_;
  print("Node types learner is using (extracted from input file):\n");
  foreach my $node (sort keys (%$ref_hash_node_types)){
    print("$node\n");
  }
}
