%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3DC*DY3DC*DZ3DC;
%Allocate space for large sparse matrix A
A = spalloc(ne,ne,7*ne);
%A = zeros(ne,ne);
    %-------------for a 3D matrix only-------------------------------------
    %----solve for global grid elements
            %----------------defines boundaries of global matrix ---------
%Activate 4 cores for parallel processing
spmd(4)
%Assign parts of sparse matrix to 1st core
if labindex==1
   for i=1:floor(DX3DC/4)
       for j = 1:DY3DC
           for k = 1:DZ3DC
%determine position within sparse matrix for each global element, using n as a pointer                
               n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;               
               %----------------defines corners of resistance matrix ---------               
               %NW corner*******from top view NW corner
               %min z
               if (i==1 && j==1) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);          
                  jpr=jparent(j,i,k);           
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));   
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  a_P  = -(a_Ex + a_Ey + a_Ez + aa_Wx + aa_Wy + aa_Wz);
                A(n,[n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_P a_Ey a_Ex a_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Ey + a_Wz + aa_Wx + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC)]) = [a_Wz a_P a_Ey a_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3DC) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wy + a_Ez + aa_Wx + aa_Ey + aa_Wz);
                A(n,[(n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wy a_P a_Ex a_Ez]; 
               %max z
               elseif (i==1 && j==DY3DC) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wy + a_Wz + aa_Wx + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC)]) = [a_Wz a_Wy a_P a_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Ex + a_Ey + a_Ez + a_Wz + aa_Wx + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_P a_Ey a_Ex a_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Ex + a_Wy + a_Ez + a_Wz + aa_Wx + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wy a_P a_Ex a_Ez];
               %W edge min z*******from top view W edge top 
               elseif (i==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Ez + aa_Wx + aa_Wz);
                A(n,[(n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wy a_P a_Ey a_Ex a_Ez];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Wz + aa_Wx + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC)]) = [a_Wz a_Wy a_P a_Ey a_Ex];
               %N edge min z******from top view N edge top
               elseif (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + aa_Wy + aa_Wz);
                A(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge bottom
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S edge top
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + aa_Ey + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               elseif i==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wx = sigma(j,i,k)*sigman(HOSTCELL-DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Ez + a_Wz + aa_Wx);
                A(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wy a_P a_Ey a_Ex a_Ez];
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];               
               end
           end
        end
   end
   for i=2:floor(DX3DC/4)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
            end
        end
    end
end
%Assign parts of sparse matrix to 2nd core
if labindex==2
   for i=floor(DX3DC/4)+1:floor(DX3DC/2)
       for j = 1:DY3DC
           for k = 1:DZ3DC
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
               %----------------defines edges of resistance matrix -----------                
               %N edge min z******from top view N edge
               if (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + aa_Wy + aa_Wz);
                A(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);         
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);        
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + aa_Ey + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/4)+1:floor(DX3DC/2)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
            end
        end
    end   
end
%Assign parts of sparse matrix to 3rd core
if labindex==3
   for i=floor(DX3DC/2)+1:floor(0.75*DX3DC)
       for j = 1:DY3DC
           for k = 1:DZ3DC
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
               %----------------defines edges of resistance matrix -----------                
               %N edge min z******from top view N edge
               if (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + aa_Wy + aa_Wz);
                A(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + aa_Ey + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/2)+1:floor(0.75*DX3DC)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
            end
        end
    end   
end
%Assign parts of sparse matrix to 4th core
if labindex==4
   for i=floor(0.75*DX3DC)+1:DX3DC
       for j = 1:DY3DC
           for k = 1:DZ3DC
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
               %NE corner******from top view NE corner
               %min z
               if (i==DX3DC && j==1) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Ey + a_Ez + aa_Ex + aa_Wy + aa_Wz);
                A(n,[(n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [a_Wx a_P a_Ey a_Ez];
               %max z
               elseif (i==DX3DC && j==1) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Ey + a_Wz + aa_Ex + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1)]) = [a_Wz a_Wx a_P a_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3DC && j==DY3DC) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);             
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Wy + a_Ez + aa_Ex + aa_Ey + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ez]; 
               %max z
               elseif (i==DX3DC && j==DY3DC) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));                
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Wy + a_Wz + aa_Ex + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n]) = [a_Wz a_Wx a_Wy a_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3DC && j==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Wx + a_Ey + a_Ez + a_Wz + aa_Ex + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_P a_Ey a_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3DC && j==DY3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Wx + a_Wy + a_Ez + a_Wz + aa_Ex + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ez];
               %E edge min z******from top view E egde  
               elseif (i==DX3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Ez + aa_Ex + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ey a_Ez];
               %E edge max z     
               elseif (i==DX3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Wz + aa_Ex + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1)]) = [a_Wz a_Wx a_Wy a_P a_Ey];
               %N edge min z******from top view N egde  
               elseif (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + aa_Wy + aa_Wz);
                A(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);            
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz + aa_Wy + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + aa_Ey + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz + aa_Ey + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               elseif i==DX3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ex = sigma(j,i,k)*sigman(HOSTCELL+DY3D)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Ez + a_Wz + aa_Ex);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ez];              
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);      
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Wy = sigma(j,i,k)*sigman(HOSTCELL-1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz + aa_Wy);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ey = sigma(j,i,k)*sigman(HOSTCELL+1)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz + aa_Ey);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  aa_Wz = sigma(j,i,k)*sigman(HOSTCELL-DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr)-zp(kpr-1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + aa_Wz);
                A(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  aa_Ez = sigma(j,i,k)*sigman(HOSTCELL+DX3D*DY3D)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(zp(kpr+2)-zp(kpr+1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz + aa_Ez);
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(0.75*DX3DC)+1:DX3DC-1
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];                                       
            end
        end
    end   
end

end
%Extract parts of sparse matrix from each core and recombine into a single
%matrix
A1 = A{:,1};
A2 = A{:,2};
A3 = A{:,3};
A4 = A{:,4};
A=A1+A2+A3+A4;