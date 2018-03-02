% AnalyzingTargetsOfDrugs.m



file = '/Users/lutzka/Dropbox/Analyses/DrugComboMethod/Data/GenericKinome/degreeinnetwork.txt';
x = importdata(file);
names = x.textdata(2:end,1);
degree = x.data;

for i = 1:length(uniqdrugs)
    for j = 1:length(sgcluster)
        substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = [];
        substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = [];
        substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = [];
        
        for k = 1:length(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).targets)
            [~,~,y] = intersect(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).targets(k),names);
            substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = ...
                [substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree; degree(y)];
        end
        for k = 1:length(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).targets)
            [~,~,y] = intersect(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).targets(k),names);
            substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = ...
                [substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree; degree(y)];
        end
        for k = 1:length(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).targets)
            [~,~,y] = intersect(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).targets(k),names);
            substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree = ...
                [substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree; degree(y)];
        end
    end
end


for i = 1:length(uniqdrugs)
    for j = 1:length(sgcluster)
        substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).maxdeg = max(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).mindeg = min(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).avgdeg = mean(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).meddeg = median(substargeted.spinglass90.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        
        substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).maxdeg = max(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).mindeg = min(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).avgdeg = mean(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).meddeg = median(substargeted.spinglass75.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        
        substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).maxdeg = max(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).mindeg = min(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).avgdeg = mean(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
        substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).meddeg = median(substargeted.spinglass50.(char(uniqdrugs(i))).(char(sgcluster(j))).degree);
    end
end



specdrugs = {};
for i = 1:(length(uniqdrugs))
    if length(inhibDrugtoConSG.spinglass50.(char(uniqdrugs(i))).clusters) == 1
        specdrugs = [specdrugs; uniqdrugs(i)];
    end
end
for i = 1:length(specdrugs)
    for j = 1:length(specdrugs)
        if i < j
            x = strcat(char(specdrugs(i)),'plus',char(specdrugs(j)));
            y = length(intersect(inhibunique50.(char(specdrugs(i))),inhibunique50.(char(specdrugs(j)))));
            z = length(union(inhibunique50.(char(specdrugs(i))),inhibunique50.(char(specdrugs(j)))));
            jacc.(char(x)) = y / z;
        end
    end
end

% inhibunique90.AZD6244 = {'MP2K1','MP2K2'};
% inhibunique75.AZD6244 = {'MP2K1','MP2K2'};
% inhibunique50.AZD6244 = {'MP2K1','MP2K2'};
% inhibunique90.JQ1 = {'BRD2','BRD3','BRD4'};
% inhibunique75.JQ1 = {'BRD2','BRD3','BRD4'};
% inhibunique50.JQ1 = {'BRD2','BRD3','BRD4'};
% uniqdrugs = [uniqdrugs; 'AZD6244';'JQ1'];
%
% x = length(intersect(inhibunique50.AZD6244,inhibunique50.Sorafenib));
% y = length(union(inhibunique50.AZD6244,inhibunique50.Sorafenib));
% jacc.azdplussoraf = x / y;
% 
% x = length(intersect(inhibunique50.AZD6244,inhibunique50.Foretinib));
% y = length(union(inhibunique50.AZD6244,inhibunique50.Foretinib));
% jacc.azdplusforet = x / y;
% 
% x = length(intersect(inhibunique50.Lapatinib,inhibunique50.JQ1));
% y = length(union(inhibunique50.Lapatinib,inhibunique50.JQ1));
% jacc.lapatplusjq1 = x / y;
% 
% x = length(intersect(inhibunique50.Sorafenib,inhibunique50.JQ1));
% y = length(union(inhibunique50.Sorafenib,inhibunique50.JQ1));
% jacc.sorafplusjq1 = x / y;
% 
% x = length(intersect(inhibunique50.Foretinib,inhibunique50.JQ1));
% y = length(union(inhibunique50.Foretinib,inhibunique50.JQ1));
% jacc.foretplusjq1 = x / y;
% 
% x = length(intersect(inhibunique50.Sorafenib,inhibunique50.Foretinib));
% y = length(union(inhibunique50.Sorafenib,inhibunique50.Foretinib));
% jacc.sorafplusforet = x / y;

for i = 1:length(uniqdrugs)
    [~,~,x] = intersect(inhibunique50.(char(uniqdrugs(i))), names);
    degreebykin50.(char(uniqdrugs(i))).mindeg = min(degree(x));
    degreebykin50.(char(uniqdrugs(i))).maxdeg = max(degree(x));
    degreebykin50.(char(uniqdrugs(i))).avgdeg = mean(degree(x));
    degreebykin50.(char(uniqdrugs(i))).meddeg = median(degree(x));
    [~,~,x] = intersect(inhibunique75.(char(uniqdrugs(i))), names);
    degreebykin75.(char(uniqdrugs(i))).mindeg = min(degree(x));
    degreebykin75.(char(uniqdrugs(i))).maxdeg = max(degree(x));
    degreebykin75.(char(uniqdrugs(i))).avgdeg = mean(degree(x));
    degreebykin75.(char(uniqdrugs(i))).meddeg = median(degree(x));
    [~,~,x] = intersect(inhibunique90.(char(uniqdrugs(i))), names);
    degreebykin90.(char(uniqdrugs(i))).mindeg = min(degree(x));
    degreebykin90.(char(uniqdrugs(i))).maxdeg = max(degree(x));
    degreebykin90.(char(uniqdrugs(i))).avgdeg = mean(degree(x));
    degreebykin90.(char(uniqdrugs(i))).meddeg = median(degree(x));
end






