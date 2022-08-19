function [] = FontWidthandPos(PosVector,FontSize,LineWidth)
if ~exist('PosVector','var') || strcmp(PosVector,'')
%     PosVector = get(gcf,'Position');
    PosVector = [50   50   800   600];
end
if ~exist('LineWidth','var') || strcmp(LineWidth,'')
    LineWidth = 2;
end
if ~exist('FontSize','var') || strcmp(FontSize,'')
    FontSize = 14;
end
%Set fontsize for all figure children
set(get(gcf,'Children'),'Fontsize',FontSize)
%Set linewidth for regular plots
set(findall(gcf,'type','line'),'LineWidth',LineWidth)
%Set linewidth for function plots
set(findall(gcf,'type','functionline'),'LineWidth',LineWidth)
%Set figure position
% set(gcf,'Units','inches','Position',PosVector)
set(gcf,'Position',PosVector)
end