%Rename a bunch of file names
paramstr = 'nu';
numold = 1;
numnew = 2;

%Get all the files that currently match
savedatapath = fullfile('F:\Dissertation\data',dirpath);
fnfrmt = @(wid,nuid,iid,raid)sprintf('data_w%s_nu%s_i%s_RA%s.mat',wid,nuid,iid,raid);

searchstr = '';
switch paramstr
    case 'w'
        searchstr = fullfile(savedatapath,fnfrmt(num2str(numold),'*','*','*'));
    case 'nu'
        searchstr = fullfile(savedatapath,fnfrmt('*',num2str(numold),'*','*'));
    case 'i'
        searchstr = fullfile(savedatapath,fnfrmt('*','*',num2str(numold),'*'));
    case 'RA'
        searchstr = fullfile(savedatapath,fnfrmt('*','*','*',num2str(numold)));
end            
fns = dir(searchstr);

for fid = 1:length(fns)
    %Make new file name
    newfn = strrep(fns(fid).name,sprintf('%s%.0f_',paramstr,numold),sprintf('%s%.0f_',paramstr,numnew));
    %Rename file
    status = movefile(fullfile(savedatapath,fns(fid).name),fullfile(savedatapath,newfn));
end