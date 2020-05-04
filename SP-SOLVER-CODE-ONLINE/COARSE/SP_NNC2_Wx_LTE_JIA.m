w = 1;
len = length(NNCC1);
while w <= len
    
    if (NNCC2(w,1)>=1  && NNCC2(w,1)<=DY3D) && (NNCC2(w,2)-NNCC1(w,2)==1)
 check = find(NNC2 == NNC2(w));
 check=check(:);
 wsize = size(check);
 nnn=1;
  for iii=1:wsize(1)
     if NNCC2(check(iii),2)-NNCC1(check(iii),2)==1
       check_cWx(nnn,1)= check(iii);
       nnn=nnn+1;
     end
  end
 
  wsize = size(check_cWx);
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
                    
            a(wcheck) = NNCC1(check_cWx(wcheck),1);
            b(wcheck) = NNCC1(check_cWx(wcheck),2);
            c(wcheck) = NNCC1(check_cWx(wcheck),3); 
%             w = w+1;            
        end
        
    else
            a(1) = NNCC1(check_cWx(1),1);
            b(1) = NNCC1(check_cWx(1),2);
            c(1) = NNCC1(check_cWx(1),3); 
            
%             w = w+1;
    end
    
%     nncnum = 0;
%     jj = 1;
%     for ii = 1:wsize(1)
% 
%         if aa(ii) == j
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
    nWx1 = (c(1)-1)*(DY3D*DX3D)+(b(1)-1)*DY3D+a(1); 
    nWx2 = (c(2)-1)*(DY3D*DX3D)+(b(2)-1)*DY3D+a(2);
    nWx2(nWx2<0)=n-DY3D;
    nWx3 = (c(3)-1)*(DY3D*DX3D)+(b(3)-1)*DY3D+a(3);
    nWx3(nWx3<0)=n-DY3D;
    nWx4 = (c(4)-1)*(DY3D*DX3D)+(b(4)-1)*DY3D+a(4);
    nWx4(nWx4<0)=n-DY3D;
    nWx5 = (c(5)-1)*(DY3D*DX3D)+(b(5)-1)*DY3D+a(5);
    nWx5(nWx5<0)=n-DY3D;
    
    D(n,(n-DY3D))=0;
%   D(n,n)=0;
    w = w+1;          
    
    
    % modifications to the loops based on the NNC 
    
    % corners 
   
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
                %NE corner******from top view NE corner
               %min z
               elseif (i==DX3D && j==1) && k==1
                  D(n,[(n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0    
                      d_Wx2 = 0;   
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(3) == 0      
                      d_Wx3 = 0;    
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(4) == 0    
                      d_Wx4 = 0;  
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(5) == 0      
                      d_Wx5 = 0;     
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Ez);
                D(n,[ n (n+1) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %max z
               elseif (i==DX3D && j==1) && k==DZ3D
                   D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1)]) = [0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  if c(2) == 0   
                      d_Wx2 = 0;   
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  end
                  if c(3) == 0     
                      d_Wx3 = 0;   
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(4) == 0        
                      d_Wx4 = 0;       
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  if c(5) == 0          
                      d_Wx5 = 0;     
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wz);
                 D(n,[(n-(DY3D*DX3D)) n (n+1)]) = [d_Wz  d_P d_Ey];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %SE corner******from top view SE corner
               %min z
               elseif (i==DX3D && j==DY3D) && k==1
                  D(n,[(n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [0 0 0 0]; 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(3) == 0       
                      d_Wx3 = 0;  
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(4) == 0     
                      d_Wx4 = 0;   
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(5) == 0       
                      d_Wx5 = 0;     
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Ez);
                D(n,[(n-1) n (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ez]; 
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5];  
               %max z
               elseif (i==DX3D && j==DY3D) && k==DZ3D
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n]) = [0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  if c(2) == 0   
                      d_Wx2 = 0;   
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(3) == 0         
                      d_Wx3 = 0;        
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(4) == 0        
                      d_Wx4 = 0;       
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  if c(5) == 0           
                      d_Wx5 = 0;          
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n]) = [d_Wz  d_Wy d_P];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
    
    % edges
    
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
               elseif (i==1 && k==1)
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez);
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wy d_P d_Ey d_Ex d_Ez];
               %W edge max z******from top view W edge bottom
               elseif (i==1 && k==DZ3D)
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D)]) = [d_Wz d_Wy d_P d_Ey d_Ex];
               %N edge min z******from top view N edge top
               elseif (j==1 && k==1)
                  D(n,[(n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0]; 
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  if c(2) == 0 
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  end
                  if c(3) == 0      
                      d_Wx3 = 0;     
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  if c(4) == 0      
                      d_Wx4 = 0;    
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(5) == 0        
                      d_Wx5 = 0;       
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Ez);
                D(n,[n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_P d_Ey d_Ex d_Ez]; 
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5];               
               %N edge max z******from top view N edge bottom
               elseif (j==1 && k==DZ3D)
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0  
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0        
                      d_Wx3 = 0;      
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  if c(4) == 0          
                      d_Wx4 = 0;    
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(5) == 0        
                      d_Wx5 = 0;       
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wz);
                 D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D)]) = [d_Wz  d_P d_Ey d_Ex];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %S edge min z******from top view S edge top
               elseif (j==DY3D && k==1)
                  D(n,[(n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0      
                      d_Wx3 = 0;     
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(4) == 0   
                      d_Wx4 = 0;  
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(5) == 0       
                      d_Wx5 = 0;      
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Ez);
                D(n,[(n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ex d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %S edge max z******from top view S edge bottom 
               elseif (j==DY3D && k==DZ3D)
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D)]) = [0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0         
                      d_Wx3 = 0;           
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);         
                  end
                  if c(4) == 0   
                      d_Wx4 = 0;   
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(5) == 0       
                      d_Wx5 = 0;      
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+DY3D)]) = [d_Wz  d_Wy d_P d_Ex];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
                %NE edge******from top view NE corner
               elseif (i==DX3D && j==1)
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  if c(2) == 0      
                      d_Wx2 = 0;     
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(3) == 0     
                      d_Wx3 = 0;    
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(4) == 0         
                      d_Wx4 = 0;         
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(5) == 0        
                      d_Wx5 = 0;         
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) n (n+1) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %SE edge******from top view SE corner 
               elseif (i==DX3D && j==DY3D)
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+(DY3D*DX3D))]) = [0 0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  if c(2) == 0    
                      d_Wx2 = 0;    
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  if c(3) == 0      
                      d_Wx3 = 0;      
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(4) == 0     
                      d_Wx4 = 0;    
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(5) == 0         
                      d_Wx5 = 0;  
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %E edge min z******from top view E egde  
               D(n,[(n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0];
               elseif (i==DX3D && k==1)
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  if c(2) == 0   
                      d_Wx2 = 0; 
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0     
                      d_Wx3 = 0;     
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(4) == 0        
                      d_Wx4 = 0;      
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  if c(5) == 0        
                      d_Wx5 = 0;      
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Ez);
                D(n,[ (n-1) n (n+1) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %E edge max z     
               elseif (i==DX3D && k==DZ3D)
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1)]) = [0 0 0 0 0];
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0   
                      d_Wx3 = 0;  
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(4) == 0     
                      d_Wx4 = 0;     
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(5) == 0         
                      d_Wx5 = 0;          
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1)]) = [d_Wz  d_Wy d_P d_Ey];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
    
    
    %faces
               %W face******from top view W face 
               elseif i==1
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz d_Wy d_P d_Ey d_Ex d_Ez];
               %N face******from top view N face
               elseif j==1
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  if c(2) == 0       
                      d_Wx2 = 0;     
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  if c(3) == 0         
                      d_Wx3 = 0;    
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(4) == 0         
                      d_Wx4 = 0;       
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  if c(5) == 0         
                      d_Wx5 = 0;        
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);          
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_P d_Ey d_Ex d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %S face******from top view S face    
               elseif j==DY3D
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  if c(2) == 0    
                      d_Wx2 = 0;   
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  if c(3) == 0      
                      d_Wx3 = 0;      
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(4) == 0          
                      d_Wx4 = 0;        
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);         
                  end
                  if c(5) == 0          
                      d_Wx5 = 0;        
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);        
                  end
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ex d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %near face******from top view top face
                elseif k==1
                  D(n,[(n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  if c(2) == 0     
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(3) == 0       
                      d_Wx3 = 0;     
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(4) == 0       
                      d_Wx4 = 0;     
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  if c(5) == 0        
                      d_Wx5 = 0;        
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Ez);
                D(n,[(n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [ d_Wy d_P d_Ey d_Ex d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
               %far face******from top view bottom face                
                elseif k==DZ3D
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D)]) = [0 0 0 0 0 0];  
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1); 
                  end
                  if c(3) == 0       
                      d_Wx3 = 0;        
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);       
                  end
                  if c(4) == 0         
                      d_Wx4 = 0;      
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  end
                  if c(5) == 0     
                      d_Wx5 = 0;   
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Wz);
                D(n,[(n-(DY3D*DX3D))  (n-1) n (n+1) (n+DY3D)]) = [d_Wz  d_Wy d_P d_Ey d_Ex];    
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5]; 
                 %E face******from top view E face
               elseif i==DX3D
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0];   
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);     
                  if c(2) == 0   
                      d_Wx2 = 0;  
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0    
                      d_Wx3 = 0;  
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(4) == 0       
                      d_Wx4 = 0; 
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(5) == 0      
                      d_Wx5 = 0;    
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);  
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Ez + d_Wz);
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ez];   
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5];   
    
    
    
    %non boundary
               else 
                  D(n,[(n-(DY3D*DX3D)) (n-DY3D) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [0 0 0 0 0 0 0];
                  d_Ex = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+2)-0.5*x(i)))*(sigma(j,i+1,k)*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(j,i+1,k)*0.5*(x(i+2)-x(i+1)))/(sigma(j,i,k)*0.5*(x(i+2)-x(i+1))+sigma(j,i+1,k)*0.5*(x(i+1)-x(i))); 
                  d_Wx1 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(1),b(1),c(1))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(1),b(1),c(1))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(1),b(1),c(1))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  if c(2) == 0   
                      d_Wx2 = 0; 
                  else
                      d_Wx2 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(2),b(2),c(2))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(2),b(2),c(2))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(2),b(2),c(2))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(3) == 0   
                      d_Wx3 = 0; 
                  else
                      d_Wx3 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(3),b(3),c(3))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(3),b(3),c(3))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(3),b(3),c(3))*0.5*(x(i+1)-x(i)))/wsize(1);   
                  end
                  if c(4) == 0        
                      d_Wx4 = 0;   
                  else
                      d_Wx4 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(4),b(4),c(4))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(4),b(4),c(4))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(4),b(4),c(4))*0.5*(x(i+1)-x(i)))/wsize(1);    
                  end
                  if c(5) == 0      
                      d_Wx5 = 0;    
                  else
                      d_Wx5 = - (((y(j+1)-y(j))*(z(k+1)-z(k)))/(0.5*x(i+1)-0.5*x(i-1)))*(sigma(a(5),b(5),c(5))*L_te(j,i,k)*0.5*(x(i+1)-x(i))+sigma(j,i,k)*L_te(a(5),b(5),c(5))*0.5*(x(i)-x(i-1)))/(sigma(j,i,k)*0.5*(x(i)-x(i-1))+sigma(a(5),b(5),c(5))*0.5*(x(i+1)-x(i)))/wsize(1);      
                  end
                  d_Ey = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+2)-0.5*y(j)))*(sigma(j+1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j+1,i,k)*0.5*(y(j+2)-y(j+1)))/(sigma(j,i,k)*0.5*(y(j+2)-y(j+1))+sigma(j+1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Wy = - (((x(i+1)-x(i))*(z(k+1)-z(k)))/(0.5*y(j+1)-0.5*y(j-1)))*(sigma(j-1,i,k)*L_te(j,i,k)*0.5*(y(j+1)-y(j))+sigma(j,i,k)*L_te(j-1,i,k)*0.5*(y(j)-y(j-1)))/(sigma(j,i,k)*0.5*(y(j)-y(j-1))+sigma(j-1,i,k)*0.5*(y(j+1)-y(j)));
                  d_Ez = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+2)-0.5*z(k)))*(sigma(j,i,k+1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k+1)*0.5*(z(k+2)-z(k+1)))/(sigma(j,i,k)*0.5*(z(k+2)-z(k+1))+sigma(j,i,k+1)*0.5*(z(k+1)-z(k)));
                  d_Wz = - (((x(i+1)-x(i))*(y(j+1)-y(j)))/(0.5*z(k+1)-0.5*z(k-1)))*(sigma(j,i,k-1)*L_te(j,i,k)*0.5*(z(k+1)-z(k))+sigma(j,i,k)*L_te(j,i,k-1)*0.5*(z(k)-z(k-1)))/(sigma(j,i,k)*0.5*(z(k)-z(k-1))+sigma(j,i,k-1)*0.5*(z(k+1)-z(k)));  
                  d_P  = -(d_Ex + d_Wx1 + d_Wx2 + d_Wx3 + d_Wx4 + d_Wx5 +  d_Ey + d_Wy + d_Ez + d_Wz);    
                D(n,[(n-(DY3D*DX3D)) (n-1) n (n+1) (n+DY3D) (n+(DY3D*DX3D))]) = [d_Wz  d_Wy d_P d_Ey d_Ex d_Ez];
                D(n,[nWx1 nWx2 nWx3 nWx4 nWx5]) = [d_Wx1 d_Wx2 d_Wx3 d_Wx4 d_Wx5];   
                
               end
    
        else w = w+1;
            
        
    end
    
end

    