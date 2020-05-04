%---------inputs-------------
% x and y grids of the well 
xsec_p=21;
ysec_p=21;

% define time steps i, state the time steps where you want cross section 
i = [48 ];

% location of the sector : a layer
zz=50;

% chang the selected sector for x and y sections
% chang the plot function for x and y axis
% chag the property you want to plot
%----------------end of inputs-------------------------

% X and Y axis. can be modified based on the grid under study
XA = zeros(length(x)-1,1);  

for cc= 1: length(x)-1
    XA(cc) = (x(cc+1)+x(cc))/2;
end

YA = zeros(length(y)-1,1);  

for cc= 1: length(y)-1
    YA(cc) = (y(cc+1)+y(cc))/2;
end



figure
% Read the selected variabled across the horizontal section selected.
% modify based on what you want to plot
for L = 1:length(i)
     LL = i(L);
     format1 = 'HOR_SEC  = Uecw%d(ysec_p,:,zz);';% can be modified based on the location of the sector required
      eval(sprintf(format1,LL));
      HOR_SEC  = HOR_SEC(:);
      plot (XA,HOR_SEC);
      title('cross section');
      hold on

     legendInfo{L} = ['TS ' num2str(LL)]; 

end


     legend(legendInfo)



 hold off

