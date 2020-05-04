%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3D*DY3D*DZ3D;
%Allocate space for large sparse matrix C
C = spalloc(ne,ne,7*ne);
%C = zeros(ne,ne);
    %-------------for a 3D matrix only-------------------------------------
    %----solve for global grid elements
            %----------------defines boundaries of global matrix ---------
%Activate 4 cores for parallel processing
spmd(4)
%Assign parts of sparse matrix to 1st core
if labindex==1
   for i=1:floor(DX3D/4)
       for j = 1:DY3D
           for k = 1:DZ3D
%determine position within sparse matrix for each global element, using n as a pointer                
               n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;               
               %----------------defines corners of resistance matrix ---------               
               %NW corner*******from top view NW corner
               %min z
               if (i==1 && j==1) && k==1
                  

                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                   
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
              
                   
                   c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3D
                 

                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
   
                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                       
                   c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D)]) = [c_Wz c_P c_Ey c_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3D) && k==1
                  

                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
 
                   
                   c_P  = -(c_Ex + c_Wy + c_Ez);
                C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez]; 
               %max z
               elseif (i==1 && j==DY3D) && k==DZ3D
                  
                   

                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   
                   c_P  = -(c_Ex + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D)]) = [c_Wz c_Wy c_P c_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   
                   c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_P c_Ey c_Ex c_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wy c_P c_Ex c_Ez];
               %W edge min z*******from top view W edge top 
              
               elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0 &&  boundaryEx(j,i)==0 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez);              
                 C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1 &&  boundaryEx(j,i)==0 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez);              
                 C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez];
                 
                  elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez);              
                 C(n,[n (n+1) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez);              
                 C(n,[(n-1) n (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez);              
                C(n,[(n-1) n (n+1) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ey c_Ez];
                
                elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==1 && boundaryEx(j,i)==0 
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez);              
                C(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ex c_Ez];
                
               elseif (i==1 && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                   
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                   
                   c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez);
                C(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ey c_Ex c_Ez];
                
                
               %W edge max z******from top view W edge bottom
               
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0 && boundaryEx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);              
                 C(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);              
                 C(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);              
                 C(n,[n (n+1) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);              
                 C(n,[(n-1) n (n-(DY3D*DX3D))]) = [c_Wy c_P c_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);              
                C(n,[(n-1) n (n+1) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ey c_Wz];
                
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);              
                C(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ex c_Wz];
                
               elseif (i==1 && k==DZ3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                 
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
             
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
             
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
          
                   
                   c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [c_Wz c_Wy c_P c_Ey c_Ex];
                
               %N edge min z******from top view N edge top
               
               elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                C(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez);
              C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z******from top view N edge bottom
               
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [c_Wx c_P c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
              C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S edge top
               
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                  C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                  C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                  C(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez);
                  C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge bottom 
               
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                  C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                  C(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                  C(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                  C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               
               elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wy c_P  c_Ex c_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey c_Ex c_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey  c_Ez];
                
                elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wy c_P  c_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wy c_P c_Ey  c_Ez];
                
                elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ex c_Ez];
                
               elseif i==1
                  
                   

                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wy c_P c_Ey c_Ex c_Ez];
               %N face******from top view N face
               
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey  c_Ez];
               
               elseif j==1
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face
               
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ez];
                
               elseif j==DY3D
                 
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
                
               %near face******from top view top face
               
               elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez);
                C(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez);
                C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez);
                C(n,[n  (n+(DY3D*DX3D))]) = [ c_P  c_Ez];
                
                elseif k==1
                  
                    
    
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face  
               
               elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                C(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
                C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wz);
                C(n,[n  (n-(DY3D*DX3D))]) = [ c_P  c_Wz]; 
                
                elseif k==DZ3D


                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];               
               end
           end
        end
   end
   for i=2:floor(DX3D/4)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                        
                  if boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                  C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey c_Ex c_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);               
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);        
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez + c_Wz);      
                C(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P  c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ex c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey  c_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey  c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P   c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ez];
                
                
                else
                
      
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                       end                        
            end
        end
    end
end
%Assign parts of sparse matrix to 2nd core
if labindex==2
   for i=floor(DX3D/4)+1:floor(DX3D/2)
       for j = 1:DY3D
           for k = 1:DZ3D
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
               %----------------defines edges of resistance matrix -----------                
               %N edge min z******from top view N edge
               if (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                C(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez);
              C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];                
               %N edge max z******from top view N edge
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [c_Wx c_P c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
              C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                  C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                  C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                  C(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez);
                  C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                  C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                  C(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                  C(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                  C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey  c_Ez];
               
               elseif j==1
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ez];
                
               elseif j==DY3D
                 
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez);
                C(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez);
                C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez);
                C(n,[n  (n+(DY3D*DX3D))]) = [ c_P  c_Ez];
                
                elseif k==1
                  
                    
    
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                C(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
                C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wz);
                C(n,[n  (n-(DY3D*DX3D))]) = [ c_P  c_Wz]; 
                
                elseif k==DZ3D


                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];                
               end
           end
       end
   end
   for i=floor(DX3D/4)+1:floor(DX3D/2)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  
                

                        if boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                  C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey c_Ex c_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);               
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);        
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez + c_Wz);      
                C(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P  c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ex c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey  c_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey  c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P   c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ez];
                
                
                else
                
      
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                       end
                                       
            end
        end
    end   
end
%Assign parts of sparse matrix to 3rd core
if labindex==3
   for i=floor(DX3D/2)+1:floor(0.75*DX3D)
       for j = 1:DY3D
           for k = 1:DZ3D
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
               %----------------defines edges of resistance matrix -----------                
               %N edge min z******from top view N edge
               if (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                C(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez);
              C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];               
               %N edge max z******from top view N edge
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [c_Wx c_P c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
              C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                  C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                  C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                  C(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez);
                  C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                  C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                  C(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                  C(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                  C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey  c_Ez];
               
               elseif j==1
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ez];
                
               elseif j==DY3D
                 
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez);
                C(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez);
                C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez);
                C(n,[n  (n+(DY3D*DX3D))]) = [ c_P  c_Ez];
                
                elseif k==1
                  
                    
    
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                C(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
                C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wz);
                C(n,[n  (n-(DY3D*DX3D))]) = [ c_P  c_Wz]; 
                
                elseif k==DZ3D


                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];               
               end
           end
       end
   end
   for i=floor(DX3D/2)+1:floor(0.75*DX3D)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  
                
           
                        if boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                  C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey c_Ex c_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);               
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);        
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez + c_Wz);      
                C(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P  c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ex c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey  c_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey  c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P   c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ez];
                
                
                else
                
      
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                       end
                                       
            end
        end
    end   
end
%Assign parts of sparse matrix to 4th core
if labindex==4
   for i=floor(0.75*DX3D)+1:DX3D
       for j = 1:DY3D
           for k = 1:DZ3D
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
               %NE corner******from top view NE corner
               %min z
               if (i==DX3D && j==1) && k==1
  
                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
 
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                   c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ez];
               %max z
               elseif (i==DX3D && j==1) && k==DZ3D

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
              
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)]) = [c_Wz c_Wx c_P c_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3D && j==DY3D) && k==1
 
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
 
                   c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ez]; 
               %max z
               elseif (i==DX3D && j==DY3D) && k==DZ3D

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n]) = [c_Wz c_Wx c_Wy c_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3D && j==1)
 
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey c_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3D && j==DY3D)

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ez];
               %E edge min z******from top view E egde  
               
               elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1) (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ey c_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez);
                C(n,[  n (n+1) (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1) (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                
               elseif (i==DX3D && k==1)

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                   c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey c_Ez];
               %E edge max z  
               
               elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1) (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ey c_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
                C(n,[  n (n+1) (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1) (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                
               elseif (i==DX3D && k==DZ3D)

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [c_Wz c_Wx c_Wy c_P c_Ey];
               %N edge min z******from top view N egde  
              elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                C(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez);
              C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];                  
               %N edge max z
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ey c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [c_Wx c_P c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ex c_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
              C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_P c_Ey c_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                  C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                  C(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez);
                  C(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                  C(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez);
                  C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
               %S edge max z 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                  C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                  C(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wy c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                  C(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                  C(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ex c_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                  C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               
               elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P  c_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey c_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))   n (n+1) (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey c_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey c_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ez];
                
               elseif i==DX3D
  
                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ez];              
               %N face******from top view N face
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P  c_Ex c_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ey  c_Ez];
               
               elseif j==1
                  
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_P c_Ey c_Ex c_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ex c_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ez];
                
               elseif j==DY3D
                 
                   
 
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez);
                C(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez);
                C(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez);
                C(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [c_Wx  c_P  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez);
                C(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez);
                C(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez);
                C(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ c_Wy c_P c_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez);
                C(n,[n  (n+(DY3D*DX3D))]) = [ c_P  c_Ez];
                
                elseif k==1
                  
                    
    
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez);
                C(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wz);
                C(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx  c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wy + c_Wz);
                C(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Wz);
                C(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [c_Wx c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wz);
                C(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ex + c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Ex c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wx + c_Wz);
                C(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [c_Wx  c_P  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wz);
                C(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Ey + c_Wy + c_Wz);
                C(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Ey  c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wy + c_Wz);
                C(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ c_Wy c_P c_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  c_P  = -(c_Wz);
                C(n,[n  (n-(DY3D*DX3D))]) = [ c_P  c_Wz]; 
                
                elseif k==DZ3D


                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Wz);
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex];               
               end
           end
       end
   end
   for i=floor(0.75*DX3D)+1:DX3D-1
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
  
                        if boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Ez + c_Wz);
                  C(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey c_Ex c_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P  c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Wy + c_Ez + c_Wz);               
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ey + c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P c_Ey c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Wy + c_Ez + c_Wz);    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ex c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ey + c_Wy + c_Ez + c_Wz);   
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);        
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey  c_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Ez + c_Wz);      
                C(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz   c_P  c_Ex c_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ex + c_Wx + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P  c_Ex c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wx + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [c_Wz c_Wx  c_P   c_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Ez + c_Wz);          
                C(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz   c_P c_Ey  c_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ey + c_Wy + c_Ez + c_Wz);             
                C(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P c_Ey  c_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Wy + c_Ez + c_Wz);       
                C(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [c_Wz  c_Wy c_P   c_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  c_P  = -(c_Ez + c_Wz);           
                C(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [c_Wz  c_P  c_Ez];
                
                
                else
                
      
                        c_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        c_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_ec(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_ec(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        c_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_ec(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_ec(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        c_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        c_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_ec(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_ec(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                c_P  = -(c_Ex + c_Wx + c_Ey + c_Wy + c_Ez + c_Wz);                    
                C(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [c_Wz c_Wx c_Wy c_P c_Ey c_Ex c_Ez];
                       end                                       
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