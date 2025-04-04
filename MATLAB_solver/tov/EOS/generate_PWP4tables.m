% Script to generate table from PWP4 routines

addpath ../
N = 500; % no pts
alleos = PWP4EOS('all');    
for k=1:length(alleos)
    fprintf('==> EOS %s\n',alleos{k});   
    [nb,e,p]=PWP4EOSTable(alleos{k},N);
end
