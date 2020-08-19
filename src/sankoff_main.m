function parsimony = sankoff_main(PhyloTree,ma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sankoff_main(PhyloTree,ma) calculates the parsimony score of PhyloTree by
% iteratively calling sankoff.m. It also displays the sequences at each
% node in PhyloTree if line 132 is uncommented.
%
% Input variables:
% PhyloTree: Phylogenetic Tree created using the neighbor join algorithm
% ma: multiple sequence alignment
%
% Output variables:
% parsimony: parsimony score of PhyloTree 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parse tree information
numLeaves = get(PhyloTree, 'NumLeaves');
branchNames = get (PhyloTree, 'BranchNames');
numBranches = get(PhyloTree, 'NumBranches');
pointers = get(PhyloTree,'Pointers');
 
%initialize individual nucleotide scores
A = [0 inf inf inf inf];
T = [inf 0, inf, inf, inf];
G = [inf, inf, 0, inf, inf];
C = [inf, inf, inf, 0, inf];
blank = [inf, inf, inf, inf, 0];

parsimony = 0;

len = size(ma(1).Sequence,2);

%% Creates tree
for j=1:len

    for i=1:numBranches
        % Pointers are the branch/leaf connectivity array
        l = pointers(i,1);
        r = pointers(i,2);
        
        % leaf node
        if (l <= numLeaves)
            if (r <= numLeaves)
                leafSeq(1) = ma(l).Sequence(j);
                leafSeq(2) = ma(r).Sequence(j);
                
                %leaf bool where 1 means it's a leaf and 0 means it's not a
                %leaf node
                leafBool(1) = 1;
                leafBool(2) = 1;
            % second node is at ancestor node
            else
                leafSeq(1) = ma(l).Sequence(j);
                internalSeq{2} = parseScores((r - numLeaves),:);
                leafBool(1) = 1;
                leafBool(2) = 0;
            end
        % first node is ancestor node
        elseif (l > numLeaves)
            if (r <= numLeaves)
                internalSeq{1} = parseScores((l - numLeaves),:);
                leafSeq(2) = ma(r).Sequence(j);
                leafBool(1) = 0;
                leafBool(2) = 1;
            % Both nodes are ancestors
            else
                internalSeq{1} = parseScores((l - numLeaves),:);
                internalSeq{2} = parseScores((r - numLeaves),:);
                leafBool(1) = 0;
                leafBool(2) = 0;
            end
        end
        
        for x=1:2
            if (leafBool(x))
                switch leafSeq(x)
                    % A, T, G, C, - in order
                    case 'A'
                        node{x} = A;
                    case 'T'
                        node{x} = T;
                    case 'G'
                        node{x} = G;
                    case 'C'
                        node{x} = C;
                    otherwise
                        node{x} = blank;
                end            
            end
        end
        
        % Performs Sankoff's algorithm on each node by determining if it is
        % an internal node or leaf node
        if (leafBool(1) && leafBool(2))
            parseScores(i,:) = sankoff(node{1}, node{2});
        elseif (leafBool(1) && leafBool(2) == 0)
            parseScores(i,:) = sankoff(node{1}, internalSeq{2});
        elseif (leafBool(1) == 0 && leafBool(2))
            parseScores(i,:) = sankoff(internalSeq{1}, node{2});
        else
            parseScores(i,:) = sankoff(internalSeq{1}, internalSeq{2});
        end
    end
    
    % Gets the minimum parsimony for each node
    minVal = min(parseScores,[],2);
    tempParsimony = minVal(numBranches);
    parsimony = parsimony + tempParsimony;

    for i=1:numBranches
        for k=1:5
            if (parseScores(i,k) == minVal(i))
                switch k
                    % A, T, G, C, -
                    case 1
                        seqs(i,j) = 'A';
                    case 2
                        seqs(i,j) = 'T';
                    case 3
                        seqs(i,j) = 'G';
                    case 4
                        seqs(i,j) = 'C';
                    otherwise
                        seqs(i,j) = '-';
                end            
            end
        end
    end
end
display(parsimony);
% Uncomment if you want to output the sequences at each node in the
% Phylogentic Tree
% display(seqs);
end