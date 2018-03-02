% ParseNetworkDatabases.m

clear all;
home;

%load in kinases master list alias data
file.masteraliases = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/KINASESmasterlist_w_Aliases.xlsx';
[~,t] = xlsread( char(file.masteraliases) );



%load in uniprot key data
%file.uniprotkey = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/uniprot-databaseHPA.tab';
%key.uniprot.raw = textread(file.uniprotkey,'%s','whitespace','\t');
load('/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/uniprot.mat');

uniprot.uniprotid = key.uniprot.raw(2:end,1);
uniprot.geneentrez = key.uniprot.raw(2:end,3);
%key.uniprot.genesyn = {};
% for i = 2:length(key.uniprot.raw)
%     x = strsplit( char( key.uniprot.raw(i) ) , '\t' );
%     key.uniprot.uniprotid(end+1) = key.uniprot.raw(2:end,1);
%     key.uniprot.geneentrez(end+1) = x(3);
%     %key.uniprot.genesyn( length(key.uniprot.geneentrez) ) = x(4);
% end

masteraliases.uniprot = t(2:end,1);
masteraliases.msgene = t(2:end,2);
masteraliases.rnagene = t(2:end,3);
masteraliases.entrezgene = t(2:end,15);

%assign each kinase in master list to a uniprot id in key.uniprot
masteraliases.uniprotid = masteraliases.uniprot;
for i = 1:length(masteraliases.uniprot)
    [~,~,i_uniprotkey] = intersect(masteraliases.entrezgene(i), uniprot.geneentrez);
    if ~isempty(i_uniprotkey)
        masteraliases.uniprotid(i) = uniprot.uniprotid(i_uniprotkey);
    else
        masteraliases.uniprotid(i) = {''};
    end
end
fprintf('loaded key files\n');



%%%load different networks
%load in hippie data
%file.hippie = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/hippie_current.txt';
%data.hippie.raw = textread(file.hippie,'%s','whitespace','\t'); %uniprot name
load('/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/hippie.mat');

hippie.node1 = {};
hippie.node2 = {};
for i = 1:length(data.hippie.raw)
    %x = strsplit( char( data.hippie.raw(i) ) , '\t' );
    y1 = strsplit( char( data.hippie.raw(i,1) ) , '_' );
    y2 = strsplit( char( data.hippie.raw(i,3) ) , '_' );
    z1 = strcmp(y1(1),masteraliases.uniprot);
    z2 = strcmp(y2(1),masteraliases.uniprot);
    if sum(z1) > 0
        if sum(z2) > 0
            hippie.node1(end+1) = masteraliases.uniprot(z1);
            hippie.node2(end+1) = masteraliases.uniprot(z2);
        end
    end
end
fprintf('loaded hippie data\n');

%load in hprd data
file.hprd = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/HPRD_ALL_BINARY_PROTEIN_PROTEIN_INTERACTIONS.xlsx';
[~,hprd.raw] = xlsread( char(file.hprd) ); %entrez gene name

hprd.node1 = {};
hprd.node2 = {};


for i = 1:length(hprd.raw(:,1))
    if strcmp(hprd.raw(i,1),'') == 0
        if strcmp(hprd.raw(i,4),'') == 0
            x1 = strcmp(hprd.raw(i,1),masteraliases.entrezgene);
            x2 = strcmp(hprd.raw(i,4),masteraliases.entrezgene);
            if sum(x1) > 0
                if sum(x2) > 0
                    hprd.node1(end+1) = masteraliases.uniprot(x1);
                    hprd.node2(end+1) = masteraliases.uniprot(x2);
                end
            end
        end
    end
end
fprintf('loaded hprd data\n');

%load phosphosite data
%file.phosphosite = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/phosphosite/Kinase_Substrate_Dataset';
%data.phosphosite.raw = textread(file.phosphosite,'%s'); %uniprot id
load('/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/phosphosite.mat');

phosphosite.node1 = {};
phosphosite.node2 = {};
for i = 2:length(data.phosphosite.raw)
    %x = strsplit( char( data.phosphosite.raw(i) ) , '\t' );
    y1 = strcmp(data.phosphosite.raw(i,2),masteraliases.uniprotid);
    y2 = strcmp(data.phosphosite.raw(i,7),masteraliases.uniprotid);
    if sum(y1) > 0
        if sum(y2) > 0
            phosphosite.node1(end+1) = masteraliases.uniprot(y1);
            phosphosite.node2(end+1) = masteraliases.uniprot(y2);
        end
    end
end
fprintf('loaded phosphosite data\n');

%load in reactome data
file.reactome = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/homo_sapiens.interactions.reactome.txt';
fid = fopen(file.reactome);
fgetl(fid);
tline = fgetl(fid);
reactome.raw = {};
reactome.node1 = {};
reactome.node2 = {};
while ischar(tline)
    x = strsplit(tline,'\t');
    if length(x) <= 6
        y1 = strsplit(char(x(1)),':');
        y2 = strsplit(char(x(2)),':');
    elseif length(x) > 6
        y1 = strsplit(char(x(1)),':');
        y2 = strsplit(char(x(3)),':');
    end
    z1 = strcmp(y1(2),masteraliases.uniprotid);
    z2 = strcmp(y2(2),masteraliases.uniprotid);
    
    if sum(z1) > 0
        if sum(z2) > 0
            if sum( strcmp(masteraliases.uniprot(z1),unique([reactome.node1,reactome.node2])) ) == 0
                %z1 not in data set yet
                reactome.node1(end+1) = masteraliases.uniprot(z1);
                reactome.node2(end+1) = masteraliases.uniprot(z2);
            elseif sum( strcmp(masteraliases.uniprot(z1),unique([reactome.node1,reactome.node2])) ) > 0
                %z1 is in data set
                if sum( strcmp(masteraliases.uniprot(z2),unique([reactome.node1,reactome.node2])) ) == 0
                    %z2 not in data set yet
                    reactome.node1(end+1) = masteraliases.uniprot(z1);
                    reactome.node2(end+1) = masteraliases.uniprot(z2);
                elseif sum( strcmp(masteraliases.uniprot(z2),unique([reactome.node1,reactome.node2])) ) > 0
                    %new node1 and new node2 are in final node1 or final node2
                    %check each pairing
                    %(rnode1,rnode2) to (finalnode1,finalnode2), done x1,y22
                    %(rnode2,rnode1) to (finalnode1,finalnode2), done x2,y12
                    %(rnode1,rnode2) to (finalnode2,finalnode1), done x3,y21
                    %(rnode2,rnode1) to (finalnode2,finalnode1)  done x4,y11
                    x1   = strcmp(masteraliases.uniprot(z1),reactome.node1);
                    y22  = strcmp(masteraliases.uniprot(z2),reactome.node2(x1));
                    
                    x2   = strcmp(masteraliases.uniprot(z2),reactome.node1);
                    y12  = strcmp(masteraliases.uniprot(z1),reactome.node2(x2));
                    
                    x3   = strcmp(masteraliases.uniprot(z1),reactome.node2);
                    y21  = strcmp(masteraliases.uniprot(z2),reactome.node1(x3));
                    
                    x4   = strcmp(masteraliases.uniprot(z2),reactome.node2);
                    y11  = strcmp(masteraliases.uniprot(z1),reactome.node1(x4));
                    if sum(y22) == 0
                        if sum(y12) == 0
                            if sum(y21) == 0
                                if sum(y11) == 0
                                    reactome.node1(end+1) = masteraliases.uniprot(z1);
                                    reactome.node2(end+1) = masteraliases.uniprot(z2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %data.reactome.raw = [data.reactome.raw; tline]; %uniprot id
    
    
    tline = fgetl(fid);
end
fclose(fid);

fprintf('loaded reactome data\n');

%load in i2d data
%file.i2d = '/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/i2d.2_9.Public.HUMAN.tab';
%data.i2d.raw = textread(file.i2d,'%s','whitespace','\t'); %uniprot id
load('/Users/lutzka/Dropbox/GomezLabShare/Kinome/NetworkData/i2d.mat');

i2d.node1 = {};
i2d.node2 = {};
for i = 2:length(data.i2d.raw)
    %x = strsplit( char( data.i2d.raw(i) ) , '\t' );
    y1 = strcmp(data.i2d.raw(i,2),masteraliases.uniprotid);
    y2 = strcmp(data.i2d.raw(i,3),masteraliases.uniprotid);
    if sum(y1) > 0
        if sum(y2) > 0
            i2d.node1(end+1) = masteraliases.uniprot(y1);
            i2d.node2(end+1) = masteraliases.uniprot(y2);
        end
    end
end
fprintf('loaded i2d data\n');


%%%compile all network information
network.node1 = {};
network.node2 = {};

[network] = compilenetworks(network,hippie);
[network] = compilenetworks(network,reactome);
[network] = compilenetworks(network,hprd);
[network] = compilenetworks(network,i2d);
[network] = compilenetworks(network,phosphosite);

network.num_nodes = length(unique([network.node1,network.node2]));
network.num_edges = length(network.node1);
network.per_kinome = round((network.num_nodes / 570)*100,2);

fid = fopen('Full_Kinome_Network_Compiled.txt','w');

fprintf(fid,'Compiled Kinome Network from the Hippie, Reactome, HPRD, I2D, and Phosphosite Databases\n');
fprintf(fid,'number of unique nodes = %i (percent of entire kinome = %f%%)\n',network.num_nodes,network.per_kinome);
fprintf(fid,'number of unique edges = %i\n',network.num_edges);

for i = 1:length(network.node1)
    fprintf(fid,'%s\t%s\n',char(network.node1(i)),char(network.node2(i)));
end

fclose(fid);










