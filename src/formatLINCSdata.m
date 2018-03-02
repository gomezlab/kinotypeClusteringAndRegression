% formatLINCSdata.m
%
% Author: Kyla Collins
% 
% Updated: 21 Nov 2016
%

cd ../Data/;
% load all of the LINCS database information from Elliot
file = 'LINCSCompDataV2.xlsx';
[kd,kt] = xlsread(file,'KinomeScanPercInhibition','','basic');
drugs.names = kt(2:end,1);

kinases.lincs_names = kt(2:end,3);
kinases.data = kd(:,3);
kinases.data(isnan(kinases.data)) = 100;
uniqdrugs = unique(drugs.names);
num_drugs = length(uniqdrugs);
uniqkinases = unique(kinases.lincs_names);
num_kinases = length(uniqkinases);

file = 'Key_LINCSnames.xlsx';
[~,t] = xlsread(file);
key.lincs = t(2:end,1);
key.kin = t(2:end,2);

%%%%%find aliases to match the network
%load in kinases master list alias data
masteraliases = 'KINASESmasterlist_w_Aliases.xlsx';
[~,t] = xlsread( char(masteraliases) ,'','','basic');

key.masteraliases.uniprot = t(2:end,1);
key.masteraliases.msgene = t(2:end,2);
key.masteraliases.rnagene = t(2:end,3);
key.masteraliases.entrezgene = t(2:end,15);



%assign each kinase in LINCS data to a uniprot name in key.masteraliases.uniprot
key.uniprot = key.kin;
for i = 1:length(key.kin)
    [~,~,i_uniprotkey] = intersect(key.kin(i), key.masteraliases.msgene);
    if ~isempty(i_uniprotkey)
        key.uniprot(i) = key.masteraliases.uniprot(i_uniprotkey);
    end
end


%90% inhibition
thresh = [50,75,10]; % 50%, 75%, and 90% inhibition 
for i = 1:length(uniqdrugs)
    inhibunique90.(char(uniqdrugs(i))) = {};
    inhib90.(char(uniqdrugs(i))) = {};
end

for i = 1:length(drugs.names)
    if kinases.data(i,1) < thresh(3)
        [~,~,ilincs] = intersect(kinases.lincs_names(i),key.lincs);
        inhib90.(char(drugs.names(i))) = [inhib90.(char(drugs.names(i))); kinases.lincs_names(i)];
        inhibunique90.(char(drugs.names(i))) = [inhibunique90.(char(drugs.names(i))); key.uniprot(ilincs)];
    end
end

for i = 1:length(uniqdrugs)
    inhibunique90.(char(uniqdrugs(i))) = unique(inhibunique90.(char(uniqdrugs(i))));
end


%read Kd data
file = 'LINCSCompDataV2.xlsx';
[kd,kt] = xlsread(file,'KinomeScanKD','','basic');
drugs.names = kt(2:end,1);

thresh2 = 1000; % 1000nM = 1uM

kinases.lincs_names = kt(2:end,3);
kinases.data = kd(:,3);
kinases.data(isnan(kinases.data)) = 1001; %1001nM
uniqdrugs = [uniqdrugs; unique(drugs.names)];
uniqdrugs = unique(uniqdrugs);
num_drugs = length(uniqdrugs);
uniqkinases = [uniqkinases; unique(kinases.lincs_names)];
uniqkinases = unique(uniqkinases);
num_kinases = length(uniqkinases);


%Kd of less than 1uM
for i = 1:length(uniqdrugs)
    if ~isfield(inhibunique90, char(uniqdrugs(i)))
        inhibunique90.(char(uniqdrugs(i))) = {};
        inhib90.(char(uniqdrugs(i))) = {};
    end
end

for i = 1:length(drugs.names)
    if kinases.data(i,1) < thresh2
        [~,~,ilincs] = intersect(kinases.lincs_names(i),key.lincs);
        inhib90.(char(drugs.names(i))) = [inhib90.(char(drugs.names(i))); kinases.lincs_names(i)];
        inhibunique90.(char(drugs.names(i))) = [inhibunique90.(char(drugs.names(i))); key.uniprot(ilincs)];
    end
end

for i = 1:length(uniqdrugs)
    inhibunique90.(char(uniqdrugs(i))) = unique(inhibunique90.(char(uniqdrugs(i))));
end

cd ../code/;
