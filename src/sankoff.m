function s = sankoff(left, right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sankoff(left, right) runs Sankoff's Algorithm at two nodes of
% phylogenetic tree
%
% Input variables:
% left: left node of branch
% right: right node of branch
%
% Output variables:
% s: parsimony score of the left node plus the score at the right node 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %A  T   G   C   -
    score = [0	3	4	9	8;
             3	0	2	4	8;
             4	2	0	4	8;
             9	4	4	0	8;
             8	8	8	8	8];

    % initialize score of output node
    s = [inf inf inf inf inf];
    
    % loop through all 5 different states and add the minimum over the 5
    % states of the left node and the minimum over the 5 states of the
    % right node
    for i = 1:5
        left_score = inf;
        for j = 1:5
            temp = score(i,j) + left(j);
            left_score = min(left_score, temp);
        end
        
        right_score = inf;
        for j = 1:5
            temp = score(i,j) + right(j);
            right_score = min(right_score, temp);
        end
        
        s(i) = left_score + right_score;
    end
    
end

