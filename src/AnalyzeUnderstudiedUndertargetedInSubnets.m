% AnalyzeUnderstudiedUndertargetedInSubnets.m
%
% Author: Kyla Collins
%
% Last Updated: 5 December 2017

clear all;
home;




%%% read in and format LINCS drug data
formatLINCSdata;

%%% read in and format KIN


%format weighted network and save
fid = fopen('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/Full_Kinome_Network_Compiled_weighted_pathways.txt');
[unformatnet,~] = textscan(fid,'%s %s %f');
format_network(unformatnet{1},unformatnet{2},unformatnet{3});
fclose(fid);

% read network that is in proper format
net = ReadNetwork('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/formatted_network_file.net');

%%% read in and format clusters found from consensus clustering with
%%% spinglass (must be in cluster > 90% of the runs to be grouped in that
%%% subnetwork)
%[d,t] = xlsread('consensusclusters_spinglass_greaterthan90percent.xls','','','basic');
fid = fopen('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/Kinase_all_clustering_membership_results.txt');
formatspec='%s %d %d %d %d %d %d';
[A,~] = textscan(fid,formatspec,'HeaderLines',1);
fclose(fid);

cd ../../code/;
n = A{1};
net.subnets.fastgreedy = [];
net.subnets.spinglass = [];
net.subnets.walktrap = [];
net.subnets.leadingeig = [];
net.subnets.labelprop = [];
net.subnets.mcode = [];
for i = 1:length(n)
    x = FindNode(net,char(n(i)));
    if x ~= -1
        %net.subnets.node(i) = net.node(i).name;
        net.subnets.fastgreedy = [net.subnets.fastgreedy; A{2}(i)];
        net.subnets.spinglass  = [net.subnets.spinglass; A{3}(i)];
        net.subnets.walktrap   = [net.subnets.walktrap; A{4}(i)];
        net.subnets.leadingeig = [net.subnets.leadingeig; A{5}(i)];
        net.subnets.labelprop  = [net.subnets.labelprop; A{6}(i)];
        net.subnets.mcode      = [net.subnets.mcode; A{7}(i)];
    end
end

%%% read in understudied kinases and format
[~,underg] = xlsread('/Users/lutzka/Dropbox/GomezLabShare/Kinome/UnderstudiedKinases.xlsx');
[~,keyu] = xlsread('/Users/lutzka/Desktop/KinomeRenderFiles/KinomeRenderTreeNames.xlsx');
keyg = keyu(2:end,1);
keyk = keyu(2:end,2);
underg = underg(2:end,1);

%convert all understudied gene names to kinase names
for i = 1:length(keyg)
    x = strcmp(keyg(i),underg);
    if sum(x) > 0
        y = find(x == 1);
        underk(y,1) = keyk(i);
    end
end

for i = 1:length(net.node)
    x = intersect(net.node(i).name,underk);
    if ~isempty(x)
        net.node(i).understudied = 1;
    elseif isempty(x)
        net.node(i).understudied = 0;
    end
end



%determine how many drugs target each individual kinase
for i = 1:length(net.node)
    net.node(i).numdrugs = 0;
end
for j = 1:length(uniqdrugs)
    for i = 1:length(net.node)
        x = intersect(net.node(i).name,inhibunique90.(char(uniqdrugs(j))));
        if ~isempty(x)
            net.node(i).numdrugs = net.node(i).numdrugs + 1;
        end
    end
end

%find degree of each node
for i = 1:length(net.node)
    net.node(i).degree = 0;
end

for j = 1:length(net.arc)
    net.node(net.arc(j).tail).degree = net.node(net.arc(j).tail).degree + 1;
    net.node(net.arc(j).head).degree = net.node(net.arc(j).head).degree + 1;
end

%analyze each clustering method individually
fastgreedyclusts = clust_stats(net, 'fastgreedy', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(fastgreedyclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/fastgreedy_stats.txt');

spinglassclusts = clust_stats(net, 'spinglass', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(spinglassclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/spinglass_stats.txt');

walktrapclusts = clust_stats(net, 'walktrap', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(walktrapclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/walktrap_stats.txt');

leadingeigclusts = clust_stats(net, 'leadingeig', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(leadingeigclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/leadingeig_stats.txt');

labelpropclusts = clust_stats(net, 'labelprop', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(labelpropclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/labelprop_stats.txt');

mcodeclusts = clust_stats(net, 'mcode', uniqdrugs, inhibunique90, n, underk);
writetable(struct2table(mcodeclusts.clusters), '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Results/mcode_stats.txt');


















