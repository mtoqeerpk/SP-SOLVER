%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3DC*DY3DC*DZ3DC;
%Allocate space for large sparse matrix B
B = spalloc(ne,ne,7*ne);
%C = zeros(ne,ne);
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
                  

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
                   
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                 
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                            bb_Wx= 0; 
                        else
                            bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                        end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                  
                   bb_Wz =0;
                   b_P  = -(b_Ex + b_Ey + b_Ez + bb_Wx + bb_Wy + bb_Wz);
                B(n,[n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_P b_Ey b_Ex b_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3DC
                 HOSTCELL = HOSTNUM_C(j,i,k);  
                 ipr=iparent(j,i,k);       
                 jpr=jparent(j,i,k);         
                 kpr=kparent(j,i,k);

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
   
                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                           bb_Wx= 0; 
                       else
                           bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                       end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k))); 
                      bb_Ez =0;
                       
                   b_P  = -(b_Ex + b_Ey + b_Wz + bb_Wx + bb_Wy + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC)]) = [b_Wz b_P b_Ey b_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3DC) && k==1
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
 
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                           bb_Wx= 0; 
                       else
                           bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                       end
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                   b_P  = -(b_Ex + b_Wy + b_Ez + bb_Wx + bb_Ey + bb_Wz);
                B(n,[(n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wy b_P b_Ex b_Ez]; 
               %max z
               elseif (i==1 && j==DY3DC) && k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
                   

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                       
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end 

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end

                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0
                           bb_Wx= 0; 
                       else
                           bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                       end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                            bb_Ey= 0; 
                        else
                            bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                        bb_Ez =0;
                   b_P  = -(b_Ex + b_Wy + b_Wz + bb_Wx  + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC)]) = [b_Wz b_Wy b_P b_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end         

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end                  

                      if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0
                          bb_Wx= 0; 
                      else
                          bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));
                      end                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                       end                 
                   
                   b_P  = -(b_Ex + b_Ey + b_Ez + b_Wz + bb_Wx + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_P b_Ey b_Ex b_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);        
                  jpr=jparent(j,i,k);      
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end                  

                      if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                          bb_Wx= 0; 
                      else
                          bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                      end
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                   b_P  = -(b_Ex + b_Wy + b_Ez + b_Wz + bb_Wx + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wy b_P b_Ex b_Ez];
               %W edge min z*******from top view W edge top 
               elseif (i==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end
                   
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
  
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                           bb_Wx= 0; 
                       else
                           bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                       end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                        b_P  = -(b_Ex + b_Ey + b_Wy + b_Ez + bb_Wx + bb_Wz);
                B(n,[(n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wy b_P b_Ey b_Ex b_Ez];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);
                  ipr=iparent(j,i,k); 
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end
                 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0;
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end
             
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end
             
                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end                  
          
                       
                      if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                          bb_Wx= 0; 
                      else
                          bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i))); 
                      end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  bb_Ez =0;
                       b_P  = -(b_Ex + b_Ey + b_Wy + b_Wz + bb_Wx + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC)]) = [b_Wz b_Wy b_P b_Ey b_Ex];
               %N edge min z******from top view N edge top
               elseif (j==1 && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end
                      
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
       
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
             
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                        bb_Wz =0;
                   
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + bb_Wy + bb_Wz);
                B(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_P b_Ey b_Ex b_Ez];               
               %N edge max z******from top view N edge bottom
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end
 
                      if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                          b_Wx = 0; 
                      else
                          b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                      end
                                          
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
                  
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                       bb_Ez =0;
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Wz + bb_Wy + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_P b_Ey b_Ex];
               %S edge min z******from top view S edge top
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   

                       if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                           b_Ex = 0; 
                       else
                           b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                       end
           
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end
            
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
                  
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                   
                   
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                            bb_Ey= 0; 
                        else
                            bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                        b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + bb_Ey + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ex b_Ez];
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   
 
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end
               
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end
            
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
              
                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
                
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                   bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Wy + b_Wz + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               elseif i==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
                  
                      if (L_ek(j,i,k)*L_ekn(HOSTCELL-DY3D))==0 
                          bb_Wx= 0; 
                      else
                          bb_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr)-xp(ipr-1))))*(sigman(HOSTCELL-DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL-DY3D)*0.5*(xp(ipr)-xp(ipr-1)))/(sigma(j,i,k)*0.5*(xp(ipr)-xp(ipr-1))+sigman(HOSTCELL-DY3D)*0.5*(x(i+1)-x(i)));
                      end
                   b_P  = -(b_Ex + b_Ey + b_Wy + b_Ez + b_Wz + bb_Wx);
                B(n,[(n-(DY3DC*DX3DC)) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wy b_P b_Ey b_Ex b_Ez];
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   
 
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end
  
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
                    
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                       end
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + b_Wz + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_P b_Ey b_Ex b_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                 HOSTCELL = HOSTNUM_C(j,i,k); 
                 ipr=iparent(j,i,k);      
                 jpr=jparent(j,i,k);   
                 kpr=kparent(j,i,k);
                   
 
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
                  
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
               
                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end
            
                         if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                             bb_Ey= 0; 
                         else
                             bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                         end
                   b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + b_Wz + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ex b_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);    
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                    
    
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
                   
                    
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                    HOSTCELL = HOSTNUM_C(j,i,k);  
                    ipr=iparent(j,i,k);  
                    jpr=jparent(j,i,k);   
                    kpr=kparent(j,i,k);

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0;
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
                
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end
                  
                        %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                    bb_Ez =0;
                        b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Wz + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex];               
               end
           end
        end
   end
   for i=2:floor(DX3DC/4)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  
                
      
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                      if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                          b_Wz = 0; 
                      else
                          b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                      end
                   
                
                b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + b_Wz);                    
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
                                       
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
                   

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                 
                   
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0;
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + bb_Wy + bb_Wz);
                B(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_P b_Ey b_Ex b_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   
          
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0;
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
                  
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wz + bb_Wy + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_P b_Ey b_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                   
      
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                    
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                        b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + bb_Ey + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ex b_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end
                  
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                 bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Wy + b_Wz + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                   
  
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0;
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                      if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                          b_Wx = 0;
                      else
                          b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                      end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
                 
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                       end
                       
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + b_Wz + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_P b_Ey b_Ex b_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);     
                  ipr=iparent(j,i,k);       
                  jpr=jparent(j,i,k);     
                  kpr=kparent(j,i,k);
                   
  
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
                   
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                       end
                   b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + b_Wz + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ex b_Ez];
               %near face******from top view top face
                elseif k==1
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                    

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0;
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
                  
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end
  
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k); 
                  ipr=iparent(j,i,k);  
                  jpr=jparent(j,i,k);  
                  kpr=kparent(j,i,k);
                    

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
                
                    
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                   bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Wz + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/4)+1:floor(DX3DC/2)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  
                

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end                 

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end         

                      if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                          b_Wz = 0; 
                      else
                          b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                      end                  
                   
                
                b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + b_Wz);                    
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
                                       
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
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end                 

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end         
                     
                   
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                        end
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + bb_Wy + bb_Wz);
                B(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_P b_Ey b_Ex b_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3DC)
                  HOSTCELL = HOSTNUM_C(j,i,k);    
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);    
                  kpr=kparent(j,i,k);

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end                 

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end                  
                   
                   
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wz + bb_Wy + bb_Ez);
                   B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_P b_Ey b_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3DC && k==1)
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k); 
                  kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                            bb_Ey= 0; 
                        else
                            bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                        b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + bb_Ey + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ex b_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3DC && k==DZ3DC)
                    HOSTCELL = HOSTNUM_C(j,i,k); 
                    ipr=iparent(j,i,k);  
                    jpr=jparent(j,i,k);
                    kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                 bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Wy + b_Wz + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1
                      HOSTCELL = HOSTNUM_C(j,i,k);   
                      ipr=iparent(j,i,k);   
                      jpr=jparent(j,i,k);  
                      kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end
                 
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                        end
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + b_Wz + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_P b_Ey b_Ex b_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);   
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);       
                  kpr=kparent(j,i,k);
                   

                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                         if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                             b_Ez = 0; 
                         else
                             b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                         end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
 
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                   b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + b_Wz + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ex b_Ez];
               %near face******from top view top face
                elseif k==1
                 
                    HOSTCELL = HOSTNUM_C(j,i,k);  
                    ipr=iparent(j,i,k);   
                    jpr=jparent(j,i,k); 
                    kpr=kparent(j,i,k);
    
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0;
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end
    
                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end
 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
                
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);   
                  jpr=jparent(j,i,k);   
                  kpr=kparent(j,i,k);
     
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0;
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0;
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
                 
                    
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                   bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Wz + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex];                
               end
           end
       end
   end
   for i=floor(DX3DC/2)+1:floor(0.75*DX3DC)
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
                  
                
           
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end
  
                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end
    
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0;
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end
   
                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end
                 
                
                b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + b_Wz);                    
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
                                       
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
                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0;
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end
 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                            bb_Ex = 0; 
                        else
                            bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i))); 
                        end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                        end
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                         b_P  = -(b_Wx + b_Ey + b_Ez + bb_Ex + bb_Wy + bb_Wz);
                B(n,[(n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [b_Wx b_P b_Ey b_Ez];
               %max z
               elseif (i==DX3DC && j==1) && k==DZ3DC
                      HOSTCELL = HOSTNUM_C(j,i,k);   
                      ipr=iparent(j,i,k);     
                      jpr=jparent(j,i,k);  
                      kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end
              
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
 
                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
 
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                            bb_Ex = 0; 
                        else
                            bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i))); 
                        end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                 bb_Ez =0;
                        b_P  = -(b_Wx + b_Ey + b_Wz + bb_Ex + bb_Wy + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1)]) = [b_Wz b_Wx b_P b_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3DC && j==DY3DC) && k==1
                       HOSTCELL = HOSTNUM_C(j,i,k);   
                       ipr=iparent(j,i,k);   
                       jpr=jparent(j,i,k); 
                       kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
 
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                            bb_Ex = 0; 
                        else
                            bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i))); 
                        end
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                            bb_Ey= 0; 
                        else
                            bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                        end
                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                         b_P  = -(b_Wx + b_Wy + b_Ez + bb_Ex + bb_Ey + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ez]; 
               %max z
               elseif (i==DX3DC && j==DY3DC) && k==DZ3DC
                       HOSTCELL = HOSTNUM_C(j,i,k);  
                       ipr=iparent(j,i,k);   
                       jpr=jparent(j,i,k); 
                       kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end
 
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                           bb_Ex = 0; 
                       else
                           bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));
                       end
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                 bb_Ez =0;
                       b_P  = -(b_Wx + b_Wy + b_Wz + bb_Ex + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n]) = [b_Wz b_Wx b_Wy b_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3DC && j==1)
                   HOSTCELL = HOSTNUM_C(j,i,k);   
                   ipr=iparent(j,i,k);      
                   jpr=jparent(j,i,k);      
                   kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0;
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                           bb_Ex = 0; 
                       else
                           bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));
                       end
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j)));
                       end
                   b_P  = -(b_Wx + b_Ey + b_Ez + b_Wz + bb_Ex + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_P b_Ey b_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3DC && j==DY3DC)
                    HOSTCELL = HOSTNUM_C(j,i,k);  
                    ipr=iparent(j,i,k);    
                    jpr=jparent(j,i,k);   
                    kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                            bb_Ex = 0; 
                        else
                            bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i))); 
                        end
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                       end
                   b_P  = -(b_Wx + b_Wy + b_Ez + b_Wz + bb_Ex + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ez];
               %E edge min z******from top view E egde  
               elseif (i==DX3DC && k==1)
                        HOSTCELL = HOSTNUM_C(j,i,k);  
                        ipr=iparent(j,i,k);  
                        jpr=jparent(j,i,k); 
                        kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
               
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                            bb_Ex = 0; 
                        else
                            bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                 bb_Wz =0;
                        b_P  = -(b_Wx + b_Ey + b_Wy + b_Ez + bb_Ex + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ey b_Ez];
               %E edge max z     
               elseif (i==DX3DC && k==DZ3DC)
                        HOSTCELL = HOSTNUM_C(j,i,k); 
                        ipr=iparent(j,i,k);  
                        jpr=jparent(j,i,k); 
                        kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0;
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                           bb_Ex = 0; 
                       else
                           bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  bb_Ez =0;
                       b_P  = -(b_Wx + b_Ey + b_Wy + b_Wz + bb_Ex + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1)]) = [b_Wz b_Wx b_Wy b_P b_Ey];
               %N edge min z******from top view N egde  
               elseif (j==1 && k==1)
                      HOSTCELL = HOSTNUM_C(j,i,k);  
                      ipr=iparent(j,i,k);    
                      jpr=jparent(j,i,k);     
                      kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
 
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0;
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end
                
                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                        b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + bb_Wy + bb_Wz);
                B(n,[(n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_P b_Ey b_Ex b_Ez];               
               %N edge max z
               elseif (j==1 && k==DZ3DC)
                    HOSTCELL = HOSTNUM_C(j,i,k);  
                    ipr=iparent(j,i,k);     
                    jpr=jparent(j,i,k);    
                    kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end
 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                           bb_Wy= 0; 
                       else
                           bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                  bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wz + bb_Wy + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_P b_Ey b_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3DC && k==1)
                    HOSTCELL = HOSTNUM_C(j,i,k);  
                    ipr=iparent(j,i,k);     
                    jpr=jparent(j,i,k);   
                    kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end
 
                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j)));
                       end
                        %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                  bb_Wz =0;
                        b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez +  bb_Ey + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ex b_Ez];
               %S edge max z 
               elseif (j==DY3DC && k==DZ3DC)
                      HOSTCELL = HOSTNUM_C(j,i,k); 
                      ipr=iparent(j,i,k);  
                      jpr=jparent(j,i,k); 
                      kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0
                           bb_Ey= 0; 
                       else
                           bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                       end
                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                 bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Wy + b_Wz + bb_Ey + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               elseif i==DX3DC
                   HOSTCELL = HOSTNUM_C(j,i,k); 
                   ipr=iparent(j,i,k);     
                   jpr=jparent(j,i,k);    
                   kpr=kparent(j,i,k);
                       if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                           b_Wx = 0; 
                       else
                           b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                       end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end
  
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                       end
 
                       if (L_ek(j,i,k)*L_ekn(HOSTCELL+DY3D))==0 
                           bb_Ex = 0; 
                       else
                           bb_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*(x(i+1)-x(i))+0.5*(xp(ipr+2)-xp(ipr+1))))*(sigman(HOSTCELL+DY3D)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ekn(HOSTCELL+DY3D)*0.5*(xp(ipr+2)-xp(ipr+1)))/(sigma(j,i,k)*0.5*(xp(ipr+2)-xp(ipr+1))+sigman(HOSTCELL+DY3D)*0.5*(x(i+1)-x(i)));
                       end
                   b_P  = -(b_Wx + b_Ey + b_Wy + b_Ez + b_Wz + bb_Ex);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ez];              
               %N face******from top view N face
               elseif j==1
                  HOSTCELL = HOSTNUM_C(j,i,k);  
                  ipr=iparent(j,i,k);      
                  jpr=jparent(j,i,k);        
                  kpr=kparent(j,i,k);
                         if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                             b_Ex = 0; 
                         else
                             b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                         end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ekn(HOSTCELL-1))==0 
                            bb_Wy= 0; 
                        else
                            bb_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr)-yp(jpr-1))))*(sigman(HOSTCELL-1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL-1)*0.5*(yp(jpr)-yp(jpr-1)))/(sigma(j,i,k)*0.5*(yp(jpr)-yp(jpr-1))+sigman(HOSTCELL-1)*0.5*(y(j+1)-y(j))); 
                        end
                   b_P  = -(b_Ex + b_Wx + b_Ey + b_Ez + b_Wz + bb_Wy);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_P b_Ey b_Ex b_Ez];
               %S face******from top view S face    
               elseif j==DY3DC
                     HOSTCELL = HOSTNUM_C(j,i,k); 
                     ipr=iparent(j,i,k);    
                     jpr=jparent(j,i,k);    
                     kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0;
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end

                        if (L_ek(j,i,k)*L_ekn(HOSTCELL+1))==0 
                            bb_Ey= 0; 
                        else
                            bb_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*(y(j+1)-y(j))+0.5*(yp(jpr+2)-yp(jpr+1))))*(sigman(HOSTCELL+1)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ekn(HOSTCELL+1)*0.5*(yp(jpr+2)-yp(jpr+1)))/(sigma(j,i,k)*0.5*(yp(jpr+2)-yp(jpr+1))+sigman(HOSTCELL+1)*0.5*(y(j+1)-y(j))); 
                        end
                   b_P  = -(b_Ex + b_Wx + b_Wy + b_Ez + b_Wz + bb_Ey);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ex b_Ez];
               %near face******from top view top face
                elseif k==1
                     HOSTCELL = HOSTNUM_C(j,i,k);  
                     ipr=iparent(j,i,k);   
                     jpr=jparent(j,i,k);  
                     kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0 
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end

                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                         %bb_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL-DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL-DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL-DX3D*DY3D)*0.5*(z(k+1)-z(k)));                      
                   bb_Wz =0;
                         b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + bb_Wz);
                B(n,[(n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wx b_Wy b_P b_Ey b_Ex b_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3DC
                      HOSTCELL = HOSTNUM_C(j,i,k);  
                      ipr=iparent(j,i,k);   
                      jpr=jparent(j,i,k);   
                      kpr=kparent(j,i,k);
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0; 
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                        end
 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0 
                            b_Ey = 0; 
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0 
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                       if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                           b_Wz = 0; 
                       else
                           b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                       end

                       %bb_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*(z(k+1)-z(k))+0.5*(zp(2)-zp(1))))*(sigman(HOSTCELL+DX3D*DY3D)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ekn(HOSTCELL+DX3D*DY3D)*0.5*(zp(2)-zp(1)))/(sigma(j,i,k)*0.5*(zp(2)-zp(1))+sigman(HOSTCELL+DX3D*DY3D)*0.5*(z(k+1)-z(k)));  
                    bb_Ez =0;
                       b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Wz + bb_Ez);
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC)]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex];                
               end
           end
       end
   end
   for i=floor(0.75*DX3DC)+1:DX3DC-1
       for j = 2:DY3DC-1
           for k = 2:DZ3DC-1
                n=(k-1)*(DY3DC*DX3DC)+(i-1)*DY3DC+j;
                %---------this deals with all the non-boundary elements
  
                        if (L_ek(j,i,k)*L_ek(j,i+1,k))==0 
                            b_Ex = 0;
                        else
                            b_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i-1,k))==0
                            b_Wx = 0; 
                        else
                            b_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ek(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ek(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                        end
 
                        if (L_ek(j,i,k)*L_ek(j+1,i,k))==0
                            b_Ey = 0;
                        else
                            b_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                        end

                        if (L_ek(j,i,k)*L_ek(j-1,i,k))==0
                            b_Wy = 0; 
                        else
                            b_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ek(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ek(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k+1))==0 
                            b_Ez = 0; 
                        else
                            b_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                        end

                        if (L_ek(j,i,k)*L_ek(j,i,k-1))==0 
                            b_Wz = 0; 
                        else
                            b_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ek(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ek(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                        end

                b_P  = -(b_Ex + b_Wx + b_Ey + b_Wy + b_Ez + b_Wz);                    
                B(n,[(n-(DY3DC*DX3DC)) (n-DY3DC) (n-1) n (n+1) (n+DY3DC) (n+(DY3DC*DX3DC))]) = [b_Wz b_Wx b_Wy b_P b_Ey b_Ex b_Ez];                                       
            end
        end
    end   
end

end
%Extract parts of sparse matrix from each core and recombine into a single
%matrix
B1 = B{:,1};
B2 = B{:,2};
B3 = B{:,3};
B4 = B{:,4};
B=B1+B2+B3+B4;