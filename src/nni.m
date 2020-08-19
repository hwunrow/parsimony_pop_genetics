function [PhyloTree, minp, count] = nni(PhyloTree, parsimony, ma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nni(PhyloTree, parsimony, ma) attempts to minimize the parsimony score
% of PhyloTree, using the nearest-neighbor interchange algorithm
%
% Input variables: 
% PhyloTree: Phylogenetic Tree created using the neighbor join algorithm
% parsimony: parsimony score of initial tree PhyloTree before swapping
% branches
% ma: multiple sequence alignment
%
% Output variables:
% PhyloTree: Phylogenetic Tree with the minimum parsimony score after 4
% iterations
% minp: minimum parsimony score
% tracker: tracker - 1 is iteration that contains the minimum parsimony
% score
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLeaves = get(PhyloTree, 'NumLeaves');
numBranches = get(PhyloTree, 'NumBranches');
% create a row vector of the leaf names in alphabetical order.
[~,OptimalOrder] = sort(get(PhyloTree,'LeafNames'));

% initialize values
Orders(1,:) = OptimalOrder;
Trees{1} = PhyloTree;
minp = parsimony;
count = 1;

% reorder the phylogenetic tree
for i=1:numBranches/2
     [TreeReordered, OptimalOrder] = reorder(PhyloTree,OptimalOrder,'Approximate',true);
     % store reordered phylo tree
     Trees{i+1} = TreeReordered; 
     Orders(i+1,:) = OptimalOrder;    
     parsimony(i+1) = sankoff_main(Trees{i+1}, ma);
    
     % check to see if the current tree has the minimum parsimony score
     if (parsimony(i+1) < minp)
         minp = parsimony(i+1);
         count = i+1;
     end
end
display(minp);
display(count);

% checks to see if the parsimony score of minimum the original score
if (count > 1)
    PhyloTree = Trees{count};
else
    disp('No tree with lower parsimony score was found');
end

% plot tree with minimum parsimony score
figure
H = plot(PhyloTree, 'Orientation', 'top');
title('Neighbor-Joining Distance Tree of SNP data of Mitochondrial DNA using Jukes-Cantor model');
ylabel('Evolutionary distance')

% create scatter plot of leafs of phylogenetic tree with grouping by population
xcoord = H.LeafDots.XData;
ycoord = H.LeafDots.YData;
pop = cell(numLeaves,1);
group = cell(numLeaves,1);
for a = 1:numLeaves
    pop{a} = H.leafNodeLabels(a).String(end-1:end);
    switch pop{a}
     case '01'
         group{a} = 'ASW';
     case '02'
         group{a} = 'CEU';
     case '03'
         group{a} = 'CHB';
     case '04'
         group{a} = 'CHD';
     case '05'
         group{a} = 'GIH';
     case '06'
         group{a} = 'JPT';
     case '07'
         group{a} = 'LWK';
     case '08' 
         group{a} = 'MEX';
     case '09'
         group{a} = 'MKK';
     case '10'
         group{a} = 'TSI';
     case '11'
         group{a} = 'YRI';
    end
end
 color = [.75 0 .75 ; % 'almost magenta' 
          0.5451 0.2706 0.0745; % 'brown'
          1.0000 0.5490 0; % 'orange'
          0.6118 00.7843 1; % 'sky blue'
          0 0 1; % 'blue'
          0 0 0; % 'black'
          1 .8 0; % 'yellow'
          0 1 1; % 'cyan'
          0 1 0; % 'green'
          1 0 1; % 'magenta'
          1 0 0	% 'red'
         ];
 marker = 'x.^v.sdp*o+'; 
 
figure
gscatter(xcoord,ycoord,group,color, marker);
title('Neighbor-Joining Distance Tree of SNP data of Mitochondrial DNA using Jukes-Cantor model');
ylabel('Evolutionary distance')
xlabel('')
ax = gca;
ax.YDir = 'reverse';
end

