%------------------inputs-------------------------------
% x and y slice of the the well in the parent model
xsec=18;
ysec=18;

% time steps
i = [0 1 2 3 4 5 6 7 8 9 10];

% chag the property you want to plot
%----------------end of inputs-------------------------

% max z fix. acn be modified based on the grid considered
ZA = zeros(length(z)-1,1);  

for cc= 1: length(z)-1
    ZA(cc) = (z(cc+1)+z(cc))/2;
end

figure

% Read the selected variabled across the vertical section selected

for L = 1:length(i)
     LL = i(L);
     format1 = 'VER_SEC  = Uecw%d(ysec,xsec,:);';
      eval(sprintf(format1,LL));
      
      VER_SEC  = VER_SEC(:);
      
      plot (VER_SEC,ZA);
      title('cross section');
      hold on

     legendInfo{L} = ['TS ' num2str(LL)]; 

end


     legend(legendInfo)



 hold off


