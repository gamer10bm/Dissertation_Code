% [~,hostname] = system('hostname');
% hostname(regexp(hostname,'[\n]')) = [];
% switch hostname
%     case 'MillerXPS' %Personal laptop
%         base_dir = 'C:\Users\Bailey\Documents\GitHub\KUbeSat_Simulations';
%     case 'MillerPC1' %Personal Desktop
%         base_dir = 'C:\Users\bmiller\OneDrive - jayhawksupportitservices\Documents\Coding\MATLAB\KUbeSat_Simulations';
%     otherwise %KU Desktop
%         base_dir = 'E:\Documents\GitHub\KUbeSat_Simulations';
% end
% addpath(fullfile(base_dir,'KUbeSat_Subfunctions'))
% addpath(fullfile(base_dir,'3D_Shape'))
addpath(genpath('Subfunctions'))

cmdsize = matlab.desktop.commandwindow.size;
    cmdline = repmat('=',1,cmdsize(1));