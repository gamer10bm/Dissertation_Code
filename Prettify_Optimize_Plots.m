fignow = gcf;

posvec = [275         735        2278         522];
FontWidthandPos(posvec)

xlims = [0 1.5e-6];
%Adjust subplots
subplot(1,4,1)
xlim(xlims)
ylim([0 10e5])
subplot(1,4,2)
ylim([0 15])
xlim(xlims)
subplot(1,4,3)
ylim([0 0.4])
xlim(xlims)
subplot(1,4,4)
ylim([0 10e4])
xlim(xlims)