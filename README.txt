To run enter:
	align (or align_sub)
in the command window to create a phylogenetic tree for all 11 populations (or 7 subpopulations).

To compute the parsimony score enter:
	sankoff_main(PhyloTree,ma) 
in the command window. This will display the parsimony score of your tree.

To run the nearest-neighbor interchange algorithm enter
	nni(PhyloTree, parsimony, ma) or (nni_sub(PhyloTree, parsimony, ma) )
in the command window where parsimony is the integer you just computed using sankoff_main. 

The intermediate files provides the subsubtree (s) and multiple sequence alignment (ma_sub) as spoken about in the paper. 

(1)	align.m 
	Uses all the SNPs on Mitochondrial DNA and applies a neighbor-joining algorithm to construct a phylogenetic tree of all the individuals.

	Input variables: none
	HAPMAP files of SNPs on Mitochondrial DNA of 11 populations
	Assumes that sequences are located in directory '/sequences'

	Output variables:none
	displays Phylogenetic Tree and scatter plot with grouping by population


(2)	sankoff_main(PhyloTree,ma) 
	Calculates the parsimony score of PhyloTree by iteratively calling sankoff.m. It also displays the sequences at each node in PhyloTree if line 132 is uncommented.

	Input variables:
	PhyloTree: Phylogenetic Tree created using the neighbor join algorithm
	ma: multiple sequence alignment

	Output variables:
	parsimony: parsimony score of PhyloTree 


(3) sankoff(left, right) 
	Runs Sankoff's Algorithm at two nodes of phylogenetic tree

	Input variables:
	left: left node of branch
	right: right node of branch

	Output variables:
	s: parsimony score of the left node plus the score at the right node 


(4)	nni(PhyloTree, parsimony, ma) 
	attempts to minimize the parsimony score of PhyloTree, using the nearest-neighbor interchange algorithm

	Input variables: 
	PhyloTree: Phylogenetic Tree created using the neighbor join algorithm
	parsimony: parsimony score of initial tree PhyloTree before swapping
	branches
	ma: multiple sequence alignment

	Output variables:
	PhyloTree: Phylogenetic Tree with the minimum parsimony score after 4
	iterations
	minp: minimum parsimony score
	tracker: tracker - 1 is iteration that contains the minimum parsimony
	score


	The following two functions are the same as align.m and nni.m. They are just altered to produce plots of the phylogenetic trees for all populations except(ASW, YRI, LWK, and MKK). 
(5)	align_sub.m 
	Uses all the SNPs on Mitochondrial DNA and applies a neighbor-joining algorithm to construct a phylogenetic tree of all the individuals.

	Input variables: none
	HAPMAP files of SNPs on Mitochondrial DNA of 11 populations
	Assumes that sequences are located in directory '/subsequences'

	Output variables:none
	displays Phylogenetic Tree and scatter plot with grouping by population


(6)	nni(PhyloTree, parsimony, ma) 
	Attempts to minimize the parsimony score of PhyloTree, using the nearest-neighbor interchange algorithm

	Input variables: 
	PhyloTree: Phylogenetic Tree created using the neighbor join algorithm
	parsimony: parsimony score of initial tree PhyloTree before swapping
	branches
	ma: multiple sequence alignment

	Output variables:
	PhyloTree: Phylogenetic Tree with the minimum parsimony score after 4
	iterations
	minp: minimum parsimony score
	tracker: tracker - 1 is iteration that contains the minimum parsimony
	score

