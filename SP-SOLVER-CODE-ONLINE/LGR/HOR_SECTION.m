
%------------------inputs-------------------------------
% x and y slice of the the well in the LGR
xsec=18;
ysec=18;

% time steps
i = [0 1 2 3 4 5 6 7 8 9 10];

% location of the sector : a layer
zz=50;

% chang the selected sector for x and y sections
% chang the plot function for x and y axis
% chag the property you want to plot
%----------------end of inputs-------------------------

% x and y 
XAc = zeros(length(xc)-1,1);  %mid of each cell

for cc= 1: length(xc)-1
    XAc(cc) = (xc(cc+1)+xc(cc))/2;
end
 
YAc = zeros(length(yc)-1,1);  %mid of each cell

for cc= 1: length(yc)-1
    YAc(cc) = (yc(cc+1)+yc(cc))/2;
end


figure
% Read the selected variabled across the horizontal section selected
for L = 1:length(i)
     LL = i(L);
     format1 = 'HOR_SEC  = Uecc%d(ysec,:,zz);';% can be modified based on the location of the sector required
      eval(sprintf(format1,LL));
      HOR_SEC  = HOR_SEC(:);
      
      plot (XAc,HOR_SEC); % modified based on x and y
      title('cross section');
      hold on

     legendInfo{L} = ['TS ' num2str(LL)]; 

end


     legend(legendInfo)



 hold off

