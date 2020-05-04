%------------------inputs-------------------------------
% x and y slice of the the well in the LGR
xsec=18;
ysec=18;

% time steps
i = [0 1 2 3 4 5 6 7 8 9 10];

% chag the property you want to plot
%----------------end of inputs-------------------------


ZAc = zeros(length(zc)-1,1);  %mid of each cell

for cc= 1: length(zc)-1
    ZAc(cc) = (zc(cc+1)+zc(cc))/2;
end

figure

% Read the selected variabled across the horizontal section selected

for L = 1:length(i)
     LL = i(L);
     format1 = 'VER_SECc  =Uecc%d(ysec,xsec,:);';
     eval(sprintf(format1,LL));
  
      
      VER_SECc  = VER_SECc(:);
      
      plot (VER_SECc,ZAc);
      title('cross section');
      hold on

     legendInfo{L} = ['TS ' num2str(LL)]; 

end


     legend(legendInfo)



 hold off


