function [PhyloTree, minp, count] = nni_sub(PhyloTree, parsimony, ma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nni_sub(PhyloTree, parsimony, ma) attempts to minimize the parsimony score
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
for i=1:3
     [TreeReordered, OptimalOrder] = reorder(PhyloTree,OptimalOrder,'Approximate',true);
     % store reordered phylo tree
     Trees{i+1} = TreeReordered;
     Orders(i+1,:) = OptimalOrder;    
     parsimony(i+1) = sankoff_main(Trees{i+1}, ma);
    plot(TreeReordered, 'Orientation', 'top');
    title(strcat('Iteration ',num2str(i),'Parsimony = ',num2str(parsimony(i+1))));
    ylabel('Evolutionary distance')
     % check to see if the current tree has the minimum parsimony score
     if (parsimony(i+1) <= minp)
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
 for i = 1:numLeaves
    pop{i} = H.leafNodeLabels(i).String(end-1:end);
    switch pop{i}
     case '01'
         group{i} = 'CEU';
     case '02'
         group{i} = 'CHB';
     case '03'
         group{i} = 'CHD';
     case '04'
         group{i} = 'GIH';
     case '05'
         group{i} = 'JPT';
     case '06' 
         group{i} = 'MEX';
     case '07'
         group{i} = 'TSI';
    end
 end
color = 'mrgbkyc';
marker = 'o+*.sdp';
figure
gscatter(xcoord,ycoord,group,color, marker);
title('Neighbor-Joining Distance Tree of SNP data of Mitochondrial DNA using Jukes-Cantor model');
ylabel('Evolutionary distance')
xlabel('')
ax = gca;
ax.YDir = 'reverse';
end