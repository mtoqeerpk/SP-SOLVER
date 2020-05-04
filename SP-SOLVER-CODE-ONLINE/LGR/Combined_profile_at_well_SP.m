%-------------INPUTS-----------------  
% x and y slice from the LGR 
ysec=18;
xsec=18;

% define time steps i, state the time steps where you want cross section 

i=[29 30 31 32 33 34 35 36  ];

%-------------INPUTS----------------- 

% define time steps i, state the time steps where you want cross section 

TIMEPROFILE_ek=zeros(DZC,length(i));
TIMEPROFILE_ec=zeros(DZC,length(i));
TIMEPROFILE_te=zeros(DZC,length(i));
TIMEPROFILE_sp=zeros(DZC,length(i));

% Read the selected variabled across the horizontal section selected
for L = 1:length(i)
     LL = i(L);
     format1 = 'TIMEPROFILE_ek(:,L) = -1000*Uekc%d(ysec,xsec,:);';
      eval(sprintf(format1,LL));
     format1 = 'TIMEPROFILE_ec(:,L) = -1000*Uecc%d(ysec,xsec,:);';
      eval(sprintf(format1,LL));
     format1 = 'TIMEPROFILE_te(:,L) = -1000*Utec%d(ysec,xsec,:);';
      eval(sprintf(format1,LL)); 
     format1 = 'TIMEPROFILE_sp(:,L) = Sc%d(ysec,xsec,:);';
     eval(sprintf(format1,LL));
  
end 



% plotting order, you need to change the Y-axis each time, add or reduce
% based on the requirement
% figure
% plot (i,TIMEPROFILE_ek(32,:));
% hold on
% plot (i,TIMEPROFILE_ec);
% plot (i,TIMEPROFILE_sp);
% title('Time Profile at the Well');
% legend('Uek','Uec','SP');
% hold off


