function [drug,name] = drug_formatting()

file = '/Users/lutzka/Dropbox/Data/Drugs_withTargetsListed.xlsx';

[~,text] = xlsread(file);

name = text(1,2:end);

for i = 1:length(name)
    
    drug.(char(name(i))) = text( 2:end, (i+1) );
    
end

for i = 1:length(name)
    
    x = ~ strcmp(drug.(char(name(i))),'');
    drug.(char(name(i))) = drug.(char(name(i)))(x);
    
end




