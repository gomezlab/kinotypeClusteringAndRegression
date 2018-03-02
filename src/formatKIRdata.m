% formatKIRdata.m
%
% Author: Kyla Collins
% 
% Date: 14 July 2014
%

clear all;
home;

% load all of the KIR database information from the overall file downloaded 
% from paper supplemental information
kirfile = 'EntireKIRDatasetFromPaper.xls';
[kd,kt] = xlsread(kirfile);
drugs.names = kt(2,2:end);
drugs.casid = kt(3,2:end);
kinases.kir_names = kt(4:end,1);
kinases.data = kd;


% "translate" KIR names to real protein and gene names in the data
[~,t] = xlsread('KIRtoUniprotNames_key.xlsx');
names.kir  = t(2:end,1);
names.gene = t(2:end,2);
names.prot = t(2:end,3);
kinases.prot_names = {};
kinases.gene_names = {};
for i = 1:length(kinases.kir_names)
    
    [~,ikey,~] = intersect(names.kir,kinases.kir_names(i));
    kinases.prot_names = [kinases.prot_names; names.prot(ikey)];
    kinases.gene_names = [kinases.gene_names; names.gene(ikey)];
    
end

% determine inhibited/activated kinases for each drug
for k = 1:length(drugs.names)
    
    inhib.(char(drugs.names(k))).kir = {};
    inhib.(char(drugs.names(k))).prot = {};
    inhib.(char(drugs.names(k))).gene = {};
    activ.(char(drugs.names(k))).kir = {};
    activ.(char(drugs.names(k))).prot = {};
    activ.(char(drugs.names(k))).gene = {};
    
    for j = 1:length(kinases.kir_names)
        
        if kinases.data(j,k) < 75
            
            inhib.(char(drugs.names(k))).kir = [inhib.(char(drugs.names(k))).kir; kinases.kir_names(j)];
            inhib.(char(drugs.names(k))).prot = [inhib.(char(drugs.names(k))).prot; kinases.prot_names(j)];
            inhib.(char(drugs.names(k))).gene = [inhib.(char(drugs.names(k))).gene; kinases.gene_names(j)];
            
        end
        
        if kinases.data(j,k) > 125
            
            activ.(char(drugs.names(k))).kir = [activ.(char(drugs.names(k))).kir; kinases.kir_names(j)];
            activ.(char(drugs.names(k))).prot = [activ.(char(drugs.names(k))).prot; kinases.prot_names(j)];
            activ.(char(drugs.names(k))).gene = [activ.(char(drugs.names(k))).gene; kinases.gene_names(j)];
            
        end
        
    end
    
end

% store everything to use in the rest of the flow algorithm
kir.kinases = kinases;
kir.drugs = drugs;
kir.targets.inhib = inhib;
kir.targets.activ = activ;

save('KIR_database.mat','kir');