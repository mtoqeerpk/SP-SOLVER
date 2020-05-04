    
% x and y slices of the well
ysec=64;
xsec=18;


% time steps
v= [29 30 31 32 33 34 35 36];    

for j = 1:length(v)
    i = v(j);
    
formatspec = 'Uek%d(Uek%d(:)==0) = Uek%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

formatspec = 'Uekw%d(:,:,:) = Uek%d(:,:,:)-Uek%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

formatspec = 'Uec%d(Uec%d(:)==0) = Uec%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

formatspec = 'Uecw%d(:,:,:) = Uec%d(:,:,:)-Uec%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

formatspec = 'Ute%d(Ute%d(:)==0) = Ute%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

formatspec = 'Utew%d(:,:,:) = Ute%d(:,:,:)-Ute%d(ysec,xsec,1);';
eval(sprintf(formatspec,i,i,i));

end 
  
     for j = 1:length(v)  
    i = v(j);
    
formatspec = 'SPw%d(:,:,:) = Uekw%d(:,:,:) + Uecw%d(:,:,:) + Utew%d(:,:,:);';
eval(sprintf(formatspec,i,i,i,i));


end

