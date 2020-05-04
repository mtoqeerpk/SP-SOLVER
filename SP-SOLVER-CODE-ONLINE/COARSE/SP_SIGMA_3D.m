%define total number of elements in combined matrix
%DX = 10;
%DY = 4;
%DZ =2;
ne = DX3D*DY3D*DZ3D;
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
   for i=1:floor(DX3D/4)
       for j = 1:DY3D
           for k = 1:DZ3D
%determine position within sparse matrix for each global element, using n as a pointer                
               n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;               
               %----------------defines corners of resistance matrix ---------               
               %NW corner*******from top view NW corner
               %min z
               if (i==1 && j==1) && k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Ez);
                A(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_P a_Ey a_Ex a_Ez];
               %max z
               elseif (i==1 && j==1) && k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D)]) = [a_Wz a_P a_Ey a_Ex];
               %SW corner*******from top view SW corner
               %min z
               elseif (i==1 && j==DY3D) && k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy + a_Ez);
                A(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wy a_P a_Ex a_Ez]; 
               %max z
               elseif (i==1 && j==DY3D) && k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D)]) = [a_Wz a_Wy a_P a_Ex];
               %----------------defines edges of resistance matrix -----------
               %NW edge*******from top view NW corner
               elseif (i==1 && j==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_P a_Ey a_Ex a_Ez];
               %SW edge*******from top view SW corner 
               elseif (i==1 && j==DY3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wy a_P a_Ex a_Ez];
               %W edge min z*******from top view W edge top 
               elseif (i==1 && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Ez);
                A(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wy a_P a_Ey a_Ex a_Ez];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [a_Wz a_Wy a_P a_Ey a_Ex];
               %N edge min z******from top view N edge top
               elseif (j==1 && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge bottom
               elseif (j==1 && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S edge top
               elseif (j==DY3D && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3D && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries
               %W face******from top view W face 
               elseif i==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wy a_P a_Ey a_Ex a_Ez];
               %N face******from top view N face
               elseif j==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];               
               end
           end
        end
   end
   for i=2:floor(DX3D/4)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
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
               if (j==1 && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min zz******from top view S edge
               elseif (j==DY3D && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face 
               elseif j==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(DX3D/4)+1:floor(DX3D/2)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
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
               if (j==1 && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z******from top view N edge
               elseif (j==1 && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S edge
               elseif (j==DY3D && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z******from top view S edge 
               elseif (j==DY3D && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries               
               %N face******from top view N face
               elseif j==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(DX3D/2)+1:floor(0.75*DX3D)
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
                                       
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
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [a_Wx a_P a_Ey a_Ez];
               %max z
               elseif (i==DX3D && j==1) && k==DZ3D
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)]) = [a_Wz a_Wx a_P a_Ey];
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3D && j==DY3D) && k==1
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ez]; 
               %max z
               elseif (i==DX3D && j==DY3D) && k==DZ3D
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));                
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n]) = [a_Wz a_Wx a_Wy a_P];
               %----------------defines edges of resistance matrix -----------
               %NE edge******from top view NE corner
               elseif (i==DX3D && j==1)
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_P a_Ey a_Ez];
               %SE edge******from top view SE corner 
               elseif (i==DX3D && j==DY3D)
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ez];
               %E edge min z******from top view E egde  
               elseif (i==DX3D && k==1)
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ey a_Ez];
               %E edge max z     
               elseif (i==DX3D && k==DZ3D)
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [a_Wz a_Wx a_Wy a_P a_Ey];
               %N edge min z******from top view N egde  
               elseif (j==1 && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_P a_Ey a_Ex a_Ez];               
               %N edge max z
               elseif (j==1 && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_P a_Ey a_Ex];
               %S edge min z******from top view S egde
               elseif (j==DY3D && k==1)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ex a_Ez];
               %S edge max z 
               elseif (j==DY3D && k==DZ3D)
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ex];
               %----------------define surface boundaries 
               %E face******from top view E face
               elseif i==DX3D
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ez];              
               %N face******from top view N face
               elseif j==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_P a_Ey a_Ex a_Ez];
               %S face******from top view S face    
               elseif j==DY3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ex a_Ez];
               %near face******from top view top face
                elseif k==1
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez);
                A(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx a_Wy a_P a_Ey a_Ex a_Ez];
               %far face******from top view bottom face                
                elseif k==DZ3D
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex];                
               end
           end
       end
   end
   for i=floor(0.75*DX3D)+1:DX3D-1
       for j = 2:DY3D-1
           for k = 2:DZ3D-1
                n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
                %---------this deals with all the non-boundary elements
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  a_Wy = -sigma(j,i,k)*sigma(j-1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j))); 
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy + a_Ez + a_Wz);                    
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx a_Wy a_P a_Ey a_Ex a_Ez];                                       
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