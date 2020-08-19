%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align_sub.m Uses all the SNPs on Mitochondrial DNA and applies a 
% neighbor-joining algorithm to construct a phylogenetic tree of all the individuals.
%
% Input variables: none
% HAPMAP files of SNPs on Mitochondrial DNA of 11 populations
% Assumes that sequences are located in directory '/sequences'
%
% Output variables:none
% displays Phylogenetic Tree and scatter plot with grouping by population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import SNP data
listing = dir('sub_sequences');
index = 1;
numCols = zeros(1,11);
pos = zeros(1,12); % pos(i):pos(i+1) is the range that contains the individuals from population i
 for i = 3:length(listing)
    delimiter = ' ';
    fid = fopen(strcat('sequences/',listing(i).name),'rt');
    tLines = fgets(fid);
    colNames = strsplit(tLines);
    numCols(i-2) = numel(strfind(tLines,delimiter)) + 1;
    pos(i-1) = numCols(i-2)-11 + pos(i-2);
    formatSpec = '%s%s%s%d%s%s%s%s%s%s';
    for j = 1:numCols(i-2)-11
        formatSpec = strcat(formatSpec,'%s');
    end
    dataArray{i-2} = textscan(fid, strcat(formatSpec,'%[^\n]'), 'Delimiter', delimiter);
    fclose(fid);
    
    for j = 1:numCols(i-2)-11
            SNP(index).Sequence = strjoin(dataArray{i-2}{j+11},'');
            SNP(index).Header = strcat(colNames{j+11},'0',num2str(i-2)); % header includes the sample name concatenated with the population identifier
            index = index + 1;
    end
 end
 
% remove identical sequences within each population
C = cellfun(@char,{SNP.Sequence},'unif',0);
[~,idx] = unique(C);
SNP_unique = SNP(idx);

 % align sequences
 ma = multialign(SNP_unique);
 % showalignment(ma);
 
 % compute pairwise distances using Jukes-Cantor
 D = seqpdist(ma,'Method','Jukes-Cantor','Alphabet', 'NT');
 
 % create phylogenetic tree using neighbor joining algorithm
 PhyloTree = seqneighjoin(D,'equivar',SNP_unique);
 figure
 H = plot(PhyloTree, 'Orientation', 'top');
 title('Neighbor-Joining Distance Tree of SNP data of Mitochondrial DNA using Jukes-Cantor model');
 ylabel('Evolutionary distance')
 
 % create scatter plot of leafs of phylogenetic tree with grouping by population
 xcoord = H.LeafDots.XData;
 ycoord = H.LeafDots.YData;
 pop = cell(length(SNP_unique),1);
 group = cell(length(SNP_unique),1);
 for i = 1:length(SNP_unique)
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
 gscatter(xcoord,ycoord,group,color,marker);
 title('Neighbor-Joining Distance Tree of SNP data of Mitochondrial DNA using Jukes-Cantor model');
 ylabel('Evolutionary distance')
 xlabel('')
 ax = gca;
 ax.YDir = 'reverse';
 
 

%% another approach to creating Phylogenetic Tree by averaging by population the distances computed with Jukes-Cantor
%  % average distances by population (ASW, CEU, CHB, etc.)
%  pos = zeros(1,12);
%  for i = 2:length(numCols)
%     pos(i) = numCols(i-1) + pos(i-1);
%  end
%  pos(12) = length(seq);
%  
% 
%  M = length(seq);
%  x = zeros(10,10);
%  for k = 1:10
%     for j = pos(k) + 1:pos(k+1)
%          for l = k:10
%              for i = pos(l+1) + 1:pos(l+2)
%                  x(l,k) = x(l,k) + D((j-1)*(M-j/2)+i-j); %D((j-1)*(M-j/2)+i-j) is distance between i-th and j-th sequences
%              end
%          end
%     end
%  end
%  
%  % find averages
%  D_avg = zeros(1,55);
%  count = 1;
%  for j = 1:10
%      for i = j:10
%          D_avg(count) = x(i,j) / (numCols(i)*numCols(j));
%          count = count + 1;
%      end
%  end
%  
%  % create phylogenetic tree
%  PhyloTree = seqneighjoin(D_avg);
