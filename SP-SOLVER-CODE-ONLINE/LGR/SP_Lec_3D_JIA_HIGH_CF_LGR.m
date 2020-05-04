%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3DC*DY3DC*DZ3DC;
%Allocate space for large sparse matrix A
C = spalloc(ne,ne,7*ne);
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
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));    
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                     
                  c_P  = -(c_Ex + c_Ey + c_Ez + cc_Wx + cc_Wy + cc_Wz);
                C(n,[n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_P c_Ey c_Ex c_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Ey + c_Wz + cc_Wx + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC)]) = [c_Wz c_P c_Ey c_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3DC) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wy + c_Ez + cc_Wx + cc_Ey + cc_Wz);
                C(n,[(n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wy c_P c_Ex c_Ez]; 
               %max z
               elseif (i==1 && j==DY3DC) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wy + c_Wz + cc_Wx + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC)]) = [c_Wz c_Wy c_P c_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz + cc_Wx + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_P c_Ey c_Ex c_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz + cc_Wx + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wy c_P c_Ex c_Ez];
               %W edge min z*******from top view W edge top 
               elseif (i==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + cc_Wx + cc_Wz);
                C(n,[(n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wy c_P c_Ey c_Ex c_Ez];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz + cc_Wx + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC)]) = [c_Wz c_Wy c_P c_Ey c_Ex];
               %N edge min z******from top view N edge top
               elseif (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + cc_Wy + cc_Wz);
                C(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z******from top view N edge bottom
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S edge top
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + cc_Ey + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               elseif i==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));     
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz + cc_Wx);
                C(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wy c_P c_Ey c_Ex c_Ez];
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];               
               end
           end
        end
   end
   for i=2:floor(DX3DC/4)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                                       
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
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + cc_Wy + cc_Wz);
                C(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);         
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);        
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + cc_Ey + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/4)+1:floor(DX3DC/2)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                                       
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
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + cc_Wy + cc_Wz);
                C(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + cc_Ey + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/2)+1:floor(0.75*DX3DC)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                                       
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
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Wx + c_Ey + c_Ez + cc_Ex + cc_Wy + cc_Wz);
                C(n,[(n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [c_Wx c_P c_Ey c_Ez];
               %max z
               elseif (i==DX3DC && j==1) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Wx + c_Ey + c_Wz + cc_Ex + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1)]) = [c_Wz c_Wx c_P c_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3DC && j==DY3DC) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);             
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Wx + c_Wy + c_Ez + cc_Ex + cc_Ey + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ez]; 
               %max z
               elseif (i==DX3DC && j==DY3DC) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));                
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Wx + c_Wy + c_Wz + cc_Ex + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n]) = [c_Wz c_Wx c_Wy c_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3DC && j==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz + cc_Ex + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_P c_Ey c_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3DC && j==DY3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz + cc_Ex + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ez];
               %E edge min z******from top view E egde  
               elseif (i==DX3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + cc_Ex + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ey c_Ez];
               %E edge max z     
               elseif (i==DX3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz + cc_Ex + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1)]) = [c_Wz c_Wx c_Wy c_P c_Ey];
               %N edge min z******from top view N egde  
               elseif (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + cc_Wy + cc_Wz);
                C(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);            
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz + cc_Wy + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + cc_Ey + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz + cc_Ey + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               elseif i==DX3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ecn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));  
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz + cc_Ex);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ez];              
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);     
                  jpr=jparent(j,i,k);      
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz + cc_Wy);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ecn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));  
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz + cc_Ey);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  cc_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                    
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + cc_Wz);
                C(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  cc_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ecn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));   
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz + cc_Ez);
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];                
               end
           end
       end
   end
   for i=floor(0.75*DX3DC)+1:DX3DC-1
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];                                       
            end
        end
    end   
end

end
%Extract parts of sparse matrix from each core and recombine into a single
%matrix
C1 = C{:,1};
C2 = C{:,2};
C3 = C{:,3};
C4 = C{:,4};
C=C1+C2+C3+C4;