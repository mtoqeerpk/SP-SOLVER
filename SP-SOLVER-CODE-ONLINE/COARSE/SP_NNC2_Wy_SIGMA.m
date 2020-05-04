w = 1;
len = length(NNCC1);
while w <= len
    
    if (NNCC2(w,1)>=1  && NNCC2(w,1)<=DX3D) && (NNCC2(w,1)-NNCC1(w,1)==1)
 check = find(NNC2 == NNC2(w));
 check=check(:);
  wsize = size(check);
  
  nnn=1;
  for iii=1:wsize(1)
     if NNCC2(check(iii),1)-NNCC1(check(iii),1)==1
       check_cWy (nnn,1)= check(iii);
       nnn=nnn+1;
     end
  end
 
  wsize = size(check_cWy);
  
   j = NNCC2(w,1);
    i = NNCC2(w,2);
    k = NNCC2(w,3);  
    
     a = zeros(10,1);
    b = zeros(10,1);
    c = zeros(10,1);
    
    aa = zeros(10,1);
    bb = zeros(10,1);
    cc = zeros(10,1);
    
    if wsize(1) > 1
        
        for wcheck = 1:wsize(1)
                    
            a(wcheck) = NNCC1(check_cWy(wcheck),1);
            b(wcheck) = NNCC1(check_cWy(wcheck),2);
            c(wcheck) = NNCC1(check_cWy(wcheck),3); 
%             w = w+1;            
        end
        
    else
            a(1) = NNCC1(check_cWy(1),1);
            b(1) = NNCC1(check_cWy(1),2);
            c(1) = NNCC1(check_cWy(1),3); 
            
%             w = w+1;
    end
    
%     nncnum = 0;
%     jj = 1;
%     for ii = 1:wsize(1)
% 
%         if bb(ii) == i
%             a(jj) = aa(ii);
%             b(jj) = bb(ii);
%             c(jj) = cc(ii);
%             jj=jj+1;
%         else
%             nncnum = nncnum+1;
%         end
% 
%     end
%     if nncnum>5
%         nncnum=5;
%     end
    
    n=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
    nWy1 = (c(1)-1)*(DY3D*DX3D)+(b(1)-1)*DY3D+a(1); 
    nWy2 = (c(2)-1)*(DY3D*DX3D)+(b(2)-1)*DY3D+a(2);
    %nWy2(nWy2<0)=n-1;
    nWy2(nWy2<0)=n-2;
    nWy3 = (c(3)-1)*(DY3D*DX3D)+(b(3)-1)*DY3D+a(3);
    %nWy3(nWy3<0)=n-1;
    nWy3(nWy3<0)=n-2;
    nWy4 = (c(4)-1)*(DY3D*DX3D)+(b(4)-1)*DY3D+a(4);
    %nWy4(nWy4<0)=n-1;
    nWy4(nWy4<0)=n-2;
    nWy5 = (c(5)-1)*(DY3D*DX3D)+(b(5)-1)*DY3D+a(5);
    %nWy5(nWy5<0)=n-1;
    nWy5(nWy5<0)=n-2;
    
    A(n,(n-1))=0;
%   A(n,n)=0;
    w = w+1;          
    
    
    % modifications to the loops based on the NNC 
    
    % corners 
   
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
                  A(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0]; 
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0   
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0      
                     a_Wy3 = 0;    
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(4) == 0         
                     a_Wy4 = 0;         
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);         
                 end
                 if c(5) == 0          
                     a_Wy5 = 0;         
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);        
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[n (n+DY3D) (n+(DY3D*DX3D))]) = [ a_P a_Ex a_Ez]; 
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5]; 
               %max z
               elseif (i==1 && j==DY3D) && k==DZ3D
                  A(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D)]) = [0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0 
                     a_Wy2 = 0;
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0       
                     a_Wy3 = 0;       
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);             
                 end
                 if c(4) == 0            
                     a_Wy4 = 0;             
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);            
                 end
                 if c(5) == 0             
                     a_Wy5 = 0;                
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);            
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+DY3D)]) = [a_Wz  a_P a_Ex];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
                %NE corner******from top view NE corner
               %min z
               elseif (i==DX3D && j==1) && k==1
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
                  A(n,[(n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [0 0 0 0]; 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0  
                     a_Wy2 = 0; 
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0        
                     a_Wy3 = 0;       
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(4) == 0    
                     a_Wy4 = 0;   
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(5) == 0      
                     a_Wy5 = 0;   
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[(n-DY3D) n (n+(DY3D*DX3D))]) = [a_Wx  a_P a_Ez]; 
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5]; 
               %max z
               elseif (i==DX3D && j==DY3D) && k==DZ3D
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n]) = [0 0 0 0];
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 if c(2) == 0   
                     a_Wy2 = 0;   
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(3) == 0    
                     a_Wy3 = 0;   
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(4) == 0          
                     a_Wy4 = 0;   
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(5) == 0        
                     a_Wy5 = 0;       
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n]) = [a_Wz a_Wx  a_P];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
    
    % edges
    
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
                  A(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 if c(2) == 0 
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(3) == 0  
                     a_Wy3 = 0;   
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(4) == 0        
                     a_Wy4 = 0;        
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(5) == 0       
                     a_Wy5 = 0;        
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);         
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz  a_P a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %W edge min z*******from top view W edge top 
               elseif (i==1 && k==1)
                  A(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 if c(2) == 0    
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(3) == 0         
                     a_Wy3 = 0;   
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(4) == 0         
                     a_Wy4 = 0;       
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);           
                 end
                 if c(5) == 0         
                     a_Wy5 = 0;       
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);          
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[ n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ a_P a_Ey a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3D)
                  A(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 if c(2) == 0 
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(3) == 0             
                     a_Wy3 = 0;     
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(4) == 0       
                     a_Wy4 = 0;       
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);           
                 end
                 if c(5) == 0       
                     a_Wy5 = 0;     
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D)]) = [a_Wz  a_P a_Ey a_Ex];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
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
                  A(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0 
                     a_Wy2 = 0;   
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(3) == 0          
                     a_Wy3 = 0;            
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(4) == 0         
                     a_Wy4 = 0;       
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(5) == 0         
                     a_Wy5 = 0;       
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[(n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx  a_P a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3D && k==DZ3D)
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0]; 
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 if c(2) == 0  
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(3) == 0        
                     a_Wy3 = 0;        
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(4) == 0        
                     a_Wy4 = 0;       
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(5) == 0 
                     a_Wy5 = 0;  
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+DY3D)]) = [a_Wz a_Wx  a_P a_Ex];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
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
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0  
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 end
                 if c(3) == 0      
                     a_Wy3 = 0;   
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(4) == 0 
                     a_Wy4 = 0; 
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(5) == 0           
                     a_Wy5 = 0;            
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);            
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+(DY3D*DX3D))]) = [a_Wz a_Wx  a_P a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %E edge min z******from top view E egde  
               elseif (i==DX3D && k==1)
                  A(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 if c(2) == 0   
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(3) == 0 
                     a_Wy3 = 0; 
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(4) == 0       
                     a_Wy4 = 0;        
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);           
                 end
                 if c(5) == 0     
                     a_Wy5 = 0;   
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [a_Wx  a_P a_Ey a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %E edge max z     
               elseif (i==DX3D && k==DZ3D)
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1); 
                 if c(2) == 0 
                     a_Wy2 = 0;
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0     
                     a_Wy3 = 0;     
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 end
                 if c(4) == 0     
                     a_Wy4 = 0;  
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);         
                 end
                 if c(5) == 0        
                     a_Wy5 = 0;      
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);          
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)]) = [a_Wz a_Wx  a_P a_Ey];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
    
    
    %faces
               %W face******from top view W face 
               elseif i==1
                  A(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 if c(2) == 0  
                     a_Wy2 = 0; 
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);          
                 end
                 if c(3) == 0       
                     a_Wy3 = 0;     
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(4) == 0           
                     a_Wy4 = 0;     
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(5) == 0       
                     a_Wy5 = 0;        
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz  a_P a_Ey a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
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
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 if c(2) == 0  
                     a_Wy2 = 0; 
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0       
                     a_Wy3 = 0; 
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(4) == 0         
                     a_Wy4 = 0;            
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(5) == 0        
                     a_Wy5 = 0;        
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);
                  
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx  a_P a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %near face******from top view top face
                elseif k==1
                  A(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];  
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0  
                     a_Wy2 = 0;  
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 end
                 if c(3) == 0     
                     a_Wy3 = 0;     
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1); 
                 end
                 if c(4) == 0      
                     a_Wy4 = 0;      
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(5) == 0           
                     a_Wy5 = 0;  
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez);
                A(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wx  a_P a_Ey a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];
               %far face******from top view bottom face                
                elseif k==DZ3D
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);  
                 if c(2) == 0   
                     a_Wy2 = 0; 
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(3) == 0   
                     a_Wy3 = 0;  
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);        
                 end
                 if c(4) == 0        
                     a_Wy4 = 0;     
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                 if c(5) == 0           
                     a_Wy5 = 0;          
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);       
                 end
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D)  n (n+1) (n+DY3D)]) = [a_Wz a_Wx  a_P a_Ey a_Ex];  
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];  
                 %E face******from top view E face
               elseif i==DX3D
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 if c(2) == 0   
                     a_Wy2 = 0;   
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);      
                 end
                 if c(3) == 0          
                     a_Wy3 = 0;          
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(4) == 0      
                     a_Wy4 = 0;    
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                 if c(5) == 0       
                     a_Wy5 = 0;     
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [a_Wz a_Wx  a_P a_Ey a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];     
    
    
    
    %non boundary
               else 
                  A(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0 0];
                  a_Ex = -sigma(j,i,k)*sigma(j,i+1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  a_Wx = -sigma(j,i,k)*sigma(j,i-1,k)*(y(j+1)-y(j))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(j,i-1,k)*0.5*(x(i+1)-x(i)));
                  a_Ey = -sigma(j,i,k)*sigma(j+1,i,k)*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                 a_Wy1 = -sigma(j,i,k)*sigma(a(1),b(1),c(1))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(1),b(1),c(1))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 if c(2) == 0   
                     a_Wy2 = 0;   
                 else
                     a_Wy2 = -sigma(j,i,k)*sigma(a(2),b(2),c(2))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(2),b(2),c(2))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(3) == 0   
                     a_Wy3 = 0;  
                 else
                     a_Wy3 = -sigma(j,i,k)*sigma(a(3),b(3),c(3))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(3),b(3),c(3))*0.5*(y(j+1)-y(j)))/wsize(1);   
                 end
                 if c(4) == 0     
                     a_Wy4 = 0;   
                 else
                     a_Wy4 = -sigma(j,i,k)*sigma(a(4),b(4),c(4))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(4),b(4),c(4))*0.5*(y(j+1)-y(j)))/wsize(1);    
                 end
                 if c(5) == 0     
                     a_Wy5 = 0;   
                 else
                     a_Wy5 = -sigma(j,i,k)*sigma(a(5),b(5),c(5))*(x(i+1)-x(i))*(z(k+1)-z(k))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(a(5),b(5),c(5))*0.5*(y(j+1)-y(j)))/wsize(1);     
                 end
                  a_Ez = -sigma(j,i,k)*sigma(j,i,k+1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  a_Wz = -sigma(j,i,k)*sigma(j,i,k-1)*(x(i+1)-x(i))*(y(j+1)-y(j))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));
                  a_P  = -(a_Ex + a_Wx + a_Ey + a_Wy1 + a_Wy2 + a_Wy3 + a_Wy4 + a_Wy5 + a_Ez + a_Wz);                   
                A(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [a_Wz a_Wx  a_P a_Ey a_Ex a_Ez];
                A(n,[nWy1 nWy2 nWy3 nWy4 nWy5]) = [a_Wy1 a_Wy2 a_Wy3 a_Wy4 a_Wy5];       
              
               end
   
    
    
    
    
    
        else w = w+1;
            
        
    end
    
end

    