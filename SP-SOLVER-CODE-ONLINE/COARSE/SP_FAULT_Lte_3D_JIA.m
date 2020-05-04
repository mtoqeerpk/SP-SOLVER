%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3D*DY3D*DZ3D;
%Allocate space for large sparse matrix C
D = spalloc(ne,ne,7*ne);
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
                  

                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                   
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
              
                   
                   d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3D
                 

                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
   
                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                       
                   d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D)]) = [d_Wz d_P d_Ey d_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3D) && k==1
                  

                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
 
                   
                   d_P  = -(d_Ex + d_Wy + d_Ez);
                D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez]; 
               %max z
               elseif (i==1 && j==DY3D) && k==DZ3D
                  
                   

                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   
                   d_P  = -(d_Ex + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D)]) = [d_Wz d_Wy d_P d_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   
                   d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_P d_Ey d_Ex d_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P d_Ex d_Ez];
               %W edge min z*******from top view W edge top 
              
               elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0 &&  boundaryEx(j,i)==0 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez);              
                 D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1 &&  boundaryEx(j,i)==0 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez);              
                 D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez];
                 
                  elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez);              
                 D(n,[n (n+1) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez);              
                 D(n,[(n-1) n (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ez];
                 
                 elseif (i==1 && k==1) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez);              
                D(n,[(n-1) n (n+1) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ey d_Ez];
                
                elseif (i==1 && k==1) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==1 && boundaryEx(j,i)==0 
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez);              
                D(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ex d_Ez];
                
               elseif (i==1 && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                   
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                   
                   d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ey d_Ex d_Ez];
                
                
               %W edge max z******from top view W edge bottom
               
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0 && boundaryEx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);              
                 D(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);              
                 D(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);              
                 D(n,[n (n+1) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);              
                 D(n,[(n-1) n (n-(DY3D*DX3D))]) = [d_Wy d_P d_Wz];
                 
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);              
                D(n,[(n-1) n (n+1) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ey d_Wz];
                
               elseif (i==1 && k==DZ3D) &&  boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);              
                D(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ex d_Wz];
                
               elseif (i==1 && k==DZ3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                 
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
             
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
             
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
          
                   
                   d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wy d_P d_Ey d_Ex];
                
               %N edge min z******from top view N edge top
               
               elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                D(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez);
              D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];               
               %N edge max z******from top view N edge bottom
               
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [d_Wx d_P d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
              D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_P d_Ey d_Ex];
               %S edge min z******from top view S edge top
               
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                  D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                  D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                  D(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez);
                  D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
               %S edge max z******from top view S edge bottom 
               
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                  D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                  D(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                  D(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                  D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               
               elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P  d_Ex d_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey  d_Ez];
                
                elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P  d_Ez];
                
                elseif i==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P d_Ey  d_Ez];
                
                elseif i==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ex d_Ez];
                
               elseif i==1
                  
                   

                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P d_Ey d_Ex d_Ez];
               %N face******from top view N face
               
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey  d_Ez];
               
               elseif j==1
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey d_Ex d_Ez];
               %S face******from top view S face
               
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ez];
                
               elseif j==DY3D
                 
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
                
               %near face******from top view top face
               
               elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez);
                D(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez);
                D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez);
                D(n,[n  (n+(DY3D*DX3D))]) = [ d_P  d_Ez];
                
                elseif k==1
                  
                    
    
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
               %far face******from top view bottom face  
               
               elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                D(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
                D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wz);
                D(n,[n  (n-(DY3D*DX3D))]) = [ d_P  d_Wz]; 
                
                elseif k==DZ3D


                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex];               
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
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                  D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey d_Ex d_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);               
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);        
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez + d_Wz);      
                D(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P  d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ex d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey  d_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey  d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P   d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ez];
                
                
                else
                
      
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);                    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
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
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                D(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez);
              D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];                
               %N edge max z******from top view N edge
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [d_Wx d_P d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
              D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_P d_Ey d_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                  D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                  D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                  D(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez);
                  D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                  D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                  D(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                  D(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                  D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey  d_Ez];
               
               elseif j==1
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey d_Ex d_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ez];
                
               elseif j==DY3D
                 
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez);
                D(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez);
                D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez);
                D(n,[n  (n+(DY3D*DX3D))]) = [ d_P  d_Ez];
                
                elseif k==1
                  
                    
    
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                D(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
                D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wz);
                D(n,[n  (n-(DY3D*DX3D))]) = [ d_P  d_Wz]; 
                
                elseif k==DZ3D


                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex];                
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
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                  D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey d_Ex d_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);               
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);        
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez + d_Wz);      
                D(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P  d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ex d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey  d_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey  d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P   d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ez];
                
                
                else
                
      
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);                    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
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
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                D(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez);
              D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];               
               %N edge max z******from top view N edge
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [d_Wx d_P d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
              D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_P d_Ey d_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                  D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                  D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                  D(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez);
                  D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                  D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                  D(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                  D(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                  D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey  d_Ez];
               
               elseif j==1
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey d_Ex d_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ez];
                
               elseif j==DY3D
                 
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez);
                D(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez);
                D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez);
                D(n,[n  (n+(DY3D*DX3D))]) = [ d_P  d_Ez];
                
                elseif k==1
                  
                    
    
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                D(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
                D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wz);
                D(n,[n  (n-(DY3D*DX3D))]) = [ d_P  d_Wz]; 
                
                elseif k==DZ3D


                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex];               
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
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                  D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey d_Ex d_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);               
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);        
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez + d_Wz);      
                D(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P  d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ex d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey  d_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey  d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P   d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ez];
                
                
                else
                
      
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);                    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
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
  
                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
 
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                   d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ez];
               %max z
               elseif (i==DX3D && j==1) && k==DZ3D

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
              
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)]) = [d_Wz d_Wx d_P d_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3D && j==DY3D) && k==1
 
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
 
                   d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ez]; 
               %max z
               elseif (i==DX3D && j==DY3D) && k==DZ3D

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n]) = [d_Wz d_Wx d_Wy d_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3D && j==1)
 
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
 
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey d_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3D && j==DY3D)

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ez];
               %E edge min z******from top view E egde  
               
               elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1) (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ey d_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez);
                D(n,[  n (n+1) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ez];
                
                 elseif (i==DX3D && k==1) && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                
               elseif (i==DX3D && k==1)

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                   d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey d_Ez];
               %E edge max z  
               
               elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1) (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ey d_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
                D(n,[  n (n+1) (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1) (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Wz];
                
                 elseif (i==DX3D && k==DZ3D) && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                
               elseif (i==DX3D && k==DZ3D)

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  

                   d_P  = -(d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [d_Wz d_Wx d_Wy d_P d_Ey];
               %N edge min z******from top view N egde  
              elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                D(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Ez]; 
                
                elseif (j==1 && k==1) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez);
              D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez]; 
              
               elseif (j==1 && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                      
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
       
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
             
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];                  
               %N edge max z
                elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D) n (n+1) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ey d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[n (n+DY3D) (n-(DY3D*DX3D))]) = [d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n-(DY3D*DX3D))]) = [d_Wx d_P d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ex d_Wz]; 
                
               elseif (j==1 && k==DZ3D) && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];    
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
              D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Wz]; 
              
               elseif (j==1 && k==DZ3D)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
 
                      d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
                                          
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_P d_Ey d_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                  D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                  D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 &&  boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez);
                  D(n,[ n (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                  D(n,[(n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Ez];
                  
                  elseif (j==DY3D && k==1) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez);
                  D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ez];
                  
               elseif (j==DY3D && k==1)
                  
                   

                       d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
           
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                   
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
               %S edge max z 
               elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                  D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                  D(n,[(n-1) n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wy d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                  D(n,[ n (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                  D(n,[(n-DY3D)  n (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ex d_Wz];
                  
                  elseif (j==DY3D && k==DZ3D) && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                  D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Wz]; 
                  
               elseif (j==DY3D && k==DZ3D)
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
               
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
            
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
              
                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               
               elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P  d_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey d_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))   n (n+1) (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey d_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ez];
                
                elseif i==DX3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ez];
                
               elseif i==DX3D
  
                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
 
                   d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ez];              
               %N face******from top view N face
               elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1 && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P  d_Ex d_Ez];
                
                elseif j==1 && boundaryEy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey  d_Ez];
               
               elseif j==1
                  
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                       d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 
  
                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                       d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                    
                   
                   d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_P d_Ey d_Ex d_Ez];
               %S face******from top view S face    
               elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ex d_Ez];
                
                elseif j==DY3D && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ez];
                
               elseif j==DY3D
                 
                   
 
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  
                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
               
                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
            
                   d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
               %near face******from top view top face
                elseif k==1 && boundaryEy(j,i)==0 && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez);
                D(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez);
                D(n,[ n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez);
                D(n,[(n-DY3D)  n  (n+(DY3D*DX3D))]) = [d_Wx  d_P  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez);
                D(n,[ n (n+1)  (n+(DY3D*DX3D))]) = [ d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1 && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez);
                D(n,[ (n-1) n  (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ez];
                
                elseif k==1 && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez);
                D(n,[n  (n+(DY3D*DX3D))]) = [ d_P  d_Ez];
                
                elseif k==1
                  
                    
    
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         
                   
                    
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez);
                D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wz);
                D(n,[ n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx  d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wy + d_Wz);
                D(n,[ (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Wz);
                D(n,[(n-DY3D)  n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1) (n+DY3D) (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-DY3D) (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [d_Wx d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wz);
                D(n,[ n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];   
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n+DY3D) (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Ex d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx + d_Wz);
                D(n,[(n-DY3D)  n  (n-(DY3D*DX3D))]) = [d_Wx  d_P  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                 D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wz);
                D(n,[ n (n+1)  (n-(DY3D*DX3D))]) = [ d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ey + d_Wy + d_Wz);
                D(n,[ (n-1) n (n+1)  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Ey  d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];    
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wy + d_Wz);
                D(n,[ (n-1) n  (n-(DY3D*DX3D))]) = [ d_Wy d_P d_Wz];
                
                elseif k==DZ3D && boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wz);
                D(n,[n  (n-(DY3D*DX3D))]) = [ d_P  d_Wz]; 
                
                elseif k==DZ3D


                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                
                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                  
                    d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex];               
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
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Ez + d_Wz);
                  D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey d_Ex d_Ez];
                  
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P  d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1 && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0         
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Wy + d_Ez + d_Wz);               
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ey + d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P d_Ey d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Wy + d_Ez + d_Wz);    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ex d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);   
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);        
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey  d_Ez];  
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0  && boundaryWx(j,i)==1 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Ez + d_Wz);      
                D(n,[(n-(DY3D*DX3D))   n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz   d_P  d_Ex d_Ez]; 
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==0 && boundaryWx(j,i)==0 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ex + d_Wx + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n  (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P  d_Ex d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==0 
                  d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wx + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D)) (n-DY3D)  n   (n+(DY3D*DX3D))]) = [d_Wz d_Wx  d_P   d_Ez];
                
                elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Ez + d_Wz);          
                D(n,[(n-(DY3D*DX3D))   n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz   d_P d_Ey  d_Ez];
                
                 elseif boundaryEy(j,i)==0  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ey + d_Wy + d_Ez + d_Wz);             
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1)  (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey  d_Ez];
                
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==0  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Wy + d_Ez + d_Wz);       
                D(n,[(n-(DY3D*DX3D))  (n-1) n   (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P   d_Ez]; 
                    
                elseif boundaryEy(j,i)==1  && boundaryWy(j,i)==1  && boundaryEx(j,i)==1  && boundaryWx(j,i)==1 
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k))); 
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k))); 
                  d_P  = -(d_Ez + d_Wz);           
                D(n,[(n-(DY3D*DX3D))  n  (n+(DY3D*DX3D))]) = [d_Wz  d_P  d_Ez];
                
                
                else
                
      
                        d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));

                        d_Wx = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(j,i-1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i-1,k)*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));                 

                        d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));

                        d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));         

                        d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));                  
                   
                
                d_P  = -(d_Ex + d_Wx + d_Ey + d_Wy + d_Ez + d_Wz);                    
                D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wx d_Wy d_P d_Ey d_Ex d_Ez];
                       end                                       
            end
        end
    end   
end

end
%Extract parts of sparse matrix from each core and recombine into a single
%matrix
D1 = D{:,1};
D2 = D{:,2};
D3 = D{:,3};
D4 = D{:,4};
D=D1+D2+D3+D4;