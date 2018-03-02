function [f] = format_network(head,tail,weight)


f.file = '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/formatted_network_file.net';
fid = fopen(f.file,'w');
str = 'NETWORKTITLE Formatted HPRD KIN';
if fid ~=-1
    fprintf(fid,'%s\r\n',str);
end

node_list = 'NODELIST';
fprintf(fid,'%s\r\n',node_list);

% making node list
nodes = unique([head,tail]);

for i = 1:length(nodes)
    fprintf(fid,'%s\r\n',char(nodes(i)));
end
fprintf(fid,'\r\n');


%add weights


edge_list = 'EDGELIST UNDIRECTED';
fprintf(fid,'%s\t%s\r\n',edge_list,'u');
for i = 1:length(head)
    fprintf(fid,'%s\t%s\t%f\r\n',char(head(i)),char(tail(i)),weight(i));
end



fclose(fid);