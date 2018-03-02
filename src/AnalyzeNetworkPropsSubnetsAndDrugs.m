


[~,t1] = xlsread('/Users/lutzka/Dropbox/GomezLabShare/Kinome/UnderstudiedKinases.xlsx');
[~,t2] = xlsread('/Users/lutzka/Dropbox/GomezLabShare/Kinome/KINASESmasterlist_w_Aliases.xlsx');

x = t1(2:end,1);
y = t2(2:end,3);
[undergene,i_under,i_alias] = intersect(x,y);
z = t2(2:end,1);
underuniprot = z(i_alias);


for i = 1:length(sgcluster)
    subs.spinglass.genericunderstud.(char(sgcluster(i))) = ...
        intersect(underuniprot,subs.spinglass.generic.(char(sgcluster(i))));
    for j = 1:length(uniqdrugs)
        substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).undertargets = ...
            intersect(underuniprot,substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).targets);
    end
end



fid = fopen('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/drugstatistics.txt','w');

fprintf(fid,'DrugName\t');
for i = 1:length(sgcluster)
    fprintf(fid,'%s#Targets\t%s#UnderstudiedTargets\t%sPercCoverage\t%sMaxDeg\t%sAvgDeg\t%sMedDeg',...
        char(sgcluster(i)),char(sgcluster(i)),char(sgcluster(i)),char(sgcluster(i)),char(sgcluster(i)),char(sgcluster(i)));
    if i ~= length(sgcluster)
        fprintf(fid,'\t');
    end
end
fprintf(fid,'\n');
for j = 1:length(uniqdrugs)
    fprintf(fid,'%s\t',char(uniqdrugs(j)));
    for i = 1:length(sgcluster)
        fprintf(fid,'%i\t%i\t%d\t%i\t%d\t%d',...
            length(substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).targets),...
            length(substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).undertargets),...
            substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).coverage,...
            substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).maxdeg,...
            substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).avgdeg,...
            substargeted.spinglass90.(char(uniqdrugs(j))).(char(sgcluster(i))).meddeg);
        if i ~= length(sgcluster)
            fprintf(fid,'\t');
        elseif i == length(sgcluster)
            fprintf(fid,'\n');
        end
    end
end

fclose(fid);





fid = fopen('/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/subnetworkstatistics.txt','w');

fprintf(fid,'Cluster\t#KinasesInCluster\t#UnderstudiedKinasesInCluster (%%Understudied)\t#Drugs\tDrugs\n');
for i = 1:length(sgcluster)
    fprintf(fid,'%s\t%i\t%i (%d%%)\t%i\t',...
        char(sgcluster(i)),...
        length(subs.spinglass.generic.(char(sgcluster(i)))),...
        length(subs.spinglass.genericunderstud.(char(sgcluster(i)))),...
        (length(subs.spinglass.genericunderstud.(char(sgcluster(i)))) / length(subs.spinglass.generic.(char(sgcluster(i)))))*100,...
        length(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs));
    if ~isempty(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs)
        for j = 1:length(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs)
            fprintf(fid,'%s',char(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs(j)));
            if j ~= length(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs)
                fprintf(fid,', ');
            elseif j == length(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs)
                fprintf(fid,'\n');
            end
        end
    elseif isempty(inhibConSGtoDrug.spinglass90.(char(sgcluster(i))).drugs)
        fprintf(fid,'\n');
    end
end

fclose(fid);








