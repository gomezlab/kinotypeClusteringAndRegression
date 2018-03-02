% WriteWeightedGraph.m


clear all;
home;

fid = fopen('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/Full_Kinome_Network_Compiled.txt');
fgetl(fid);
fgetl(fid);
fgetl(fid);
tline = fgetl(fid);

network.node1 = {};
network.node2 = {};
while ischar(tline)
    x = strsplit(tline,'\t');
    network.node1(end+1) = x(1);
    network.node2(end+1) = x(2);
    tline = fgetl(fid);
end

fclose(fid);


[d,t] = xlsread('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/KeggPathwayScores.xlsx','Pathway');
rowWname = t(3:end,1);
colWname = t(1,3:end);
Wmat = d;

network.arcweight = [];
for i = 1:length(network.node1)
    fprintf('i=%i\n',i);
    [~,~,i1w] = intersect(network.node1(i),rowWname);
    [~,~,j2w] = intersect(network.node2(i),colWname);
    if ~isempty(i1w) && ~isempty(j2w)
        if Wmat(i1w,j2w) ~= 0
            network.arcweight(i) = Wmat(i1w,j2w);
        else
            network.arcweight(i) = 0.5;
        end
    else
        network.arcweight(i) = 0.5;
    end
end




file = '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/Full_Kinome_Network_Compiled_weighted_pathways.txt';
fid = fopen(file,'w');

for i = 1:length(network.node1)
    fprintf(fid,'%s\t%s\t%f\n',char(network.node1(i)),char(network.node2(i)),abs(network.arcweight(i)));
end

fclose(fid);