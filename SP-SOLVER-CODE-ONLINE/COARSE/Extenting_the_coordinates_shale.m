% original file from amadi is ZNNC scribt

%DX=40;
%DY=20;
%DZ=10;
%DX3D=40;
%DY3D=20;
%DZ3D=10;

%z2=zeros(6,8,4);% to be deleted 
%z2(:,:,1)=E1;z2(:,:,2)=E2;z2(:,:,3)=E3;z2(:,:,4)=E4; % Es from the excel sheet

%ignore any shale in x and y as no nnc with pass through them

Z = z2(:,:,1:2:DZ*2);% top z cordinate of each grid


%------------adjust x and y shales---------------------
Zt=zeros(2*DY3D,2*DX3D,DZ);

west_shale=0;
east_shale=0;
north_shale=5;
south_shale=0;
top_shale=20;
bottom_shale=20;


Zt(2*north_shale+1:2*(north_shale+DY),2*west_shale+1:2*(west_shale+DX),:)=Z;

Zbottom=ones(DY3D*2,DX3D*2,1);
Zbottom(2*north_shale+1:2*(north_shale+DY),2*west_shale+1:2*(west_shale+DX),1)=z2(:,:,2*DZ);

%-----------------------------------------
if north_shale==0 && south_shale==0 
    
    Zt(1,(2*west_shale+1):2*(west_shale+DX),:)=Z(1,:,:);
    Zt(2*(north_shale+DY),(2*west_shale+1):2*(west_shale+DX),:)=Z(2*DY,:,:);
    
    Zbottom(1,(2*west_shale+1):2*(west_shale+DX),1)=z2(1,:,2*DZ);
    Zbottom(2*(north_shale+DY),(2*west_shale+1):2*(west_shale+DX),:)=z2(2*DY,:,2*DZ);
    
elseif north_shale==0
    Zt(1,(2*west_shale+1):2*(west_shale+DX),:)=Z(1,:,:);
    Zbottom(1,(2*west_shale+1):2*(west_shale+DX),1)=z2(1,:,2*DZ);
    
    for i=1:2*south_shale
       Zt(2*(north_shale+DY)+i,(2*west_shale+1):2*(west_shale+DX),:)=Z(2*DY,:,:);
       Zbottom(2*(north_shale+DY)+i,(2*west_shale+1):2*(west_shale+DX),1)=z2(2*DY,:,2*DZ);
    end

elseif south_shale==0

    for i=1:2*north_shale
       Zt(i,(2*west_shale+1):2*(west_shale+DX),:)=Z(1,:,:);
       Zbottom(i,(2*west_shale+1):2*(west_shale+DX),1)=z2(1,:,2*DZ);
    end
    Zt(2*(north_shale+DY),(2*west_shale+1):2*(west_shale+DX),:)=Z(2*DY,:,:);
    Zbottom(2*(north_shale+DY),(2*west_shale+1):2*(west_shale+DX),1)=z2(2*DY,:,2*DZ);

else

    for i=1:2*south_shale
       Zt(2*(north_shale+DY)+i,(2*west_shale+1):2*(west_shale+DX),:)=Z(2*DY,:,:);
       Zbottom(2*(north_shale+DY)+i,(2*west_shale+1):2*(west_shale+DX),1)=z2(2*DY,:,2*DZ);
    end
    for i=1:2*north_shale
       Zt(i,(2*west_shale+1):2*(west_shale+DX),:)=Z(1,:,:);
       Zbottom(i,(2*west_shale+1):2*(west_shale+DX),1)=z2(1,:,2*DZ);
    end

end 
%----------------------------------------------
if west_shale==0 && east_shale==0 
    
    Zt((2*north_shale+1):2*(north_shale+DY),1,:)=Z(:,1,:);
    Zt((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX),:)=Z(:,2*DX,:);
    
    Zbottom((2*north_shale+1):2*(north_shale+DY),1,1)=z2(:,1,2*DZ);
    Zbottom((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX),1)=z2(:,2*DX,2*DZ);

elseif west_shale==0
    Zt((2*north_shale+1):2*(north_shale+DY),1,:)=Z(:,1,:);
    Zbottom((2*north_shale+1):2*(north_shale+DY),1,1)=z2(:,1,2*DZ);

    for i=1:2*east_shale
       Zt((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX)+i,:)=Z(:,2*DX,:);
       Zbottom((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX)+i,1)=z2(:,2*DX,2*DZ);
    end

elseif east_shale==0

    for i=1:2*west_shale
       Zt((2*north_shale+1):2*(north_shale+DY),i,:)=Z(:,1,:);
       Zbottom((2*north_shale+1):2*(north_shale+DY),i,1)=z2(:,1,2*DZ);
    end
    Zt((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX),:)=Z(:,2*DX,:);
    Zbottom((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX),1)=z2(:,2*DX,2*DZ);

else

    for i=1:2*east_shale
       Zt((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX)+i,:)=Z(:,2*DX,:);
       Zbottom((2*north_shale+1):2*(north_shale+DY),2*(west_shale+DX)+i,1)=z2(:,2*DX,2*DZ);
    end
    for i=1:2*west_shale
       Zt((2*north_shale+1):2*(north_shale+DY),i,:)=Z(:,1,:);
       Zbottom((2*north_shale+1):2*(north_shale+DY),i,1)=z2(:,1,2*DZ);
    end

end 
%-----------------------------------------------------

%----------adjust the corners-----------------------

if north_shale==0 && south_shale==0
    
elseif south_shale==0
    if west_shale==0
    else
        for j=1:2*north_shale
        for i=1:2*west_shale
            Zt(j,i,:)=Z(1,1,:);
            Zbottom(j,i,1)=z2(1,1,2*DZ);
        end
        end
    end
    if east_shale==0
    else 
        for j=1:2*north_shale
        for i=1:2*east_shale
            Zt(j,2*(west_shale+DX)+i,:)=Z(1,2*DX,:);
            Zbottom(j,2*(west_shale+DX)+i,1)=z2(1,2*DX,2*DZ);
        end
        end
    end
    
elseif north_shale==0
    if west_shale==0
    else
        for j=1:2*south_shale
        for i=1:2*west_shale
            Zt(2*(north_shale+DY)+j,i,:)=Z(2*DY,1,:);
            Zbottom(2*(north_shale+DY)+j,i,1)=z2(2*DY,1,2*DZ);
        end
        end
    end
    if east_shale==0
    else 
        for j=1:2*south_shale 
        for i=1:2*east_shale
            Zt(2*(north_shale+DY)+j,2*(west_shale+DX)+i,:)=Z(2*DY,2*DX,:);
            Zbottom(2*(north_shale+DY)+j,2*(west_shale+DX)+i,1)=z2(2*DY,2*DX,2*DZ);
        end
        end
    end
    
else 
    if west_shale==0
    else 
        for j=1:2*north_shale
        for i=1:2*west_shale
            Zt(j,i,:)=Z(1,1,:);
            Zbottom(j,i,1)=z2(1,1,2*DZ);
        end
        end
        for j=1:2*south_shale
        for i=1:2*west_shale
            Zt(2*(north_shale+DY)+j,i,:)=Z(2*DY,1,:);
            Zbottom(2*(north_shale+DY)+j,i,1)=z2(2*DY,1,2*DZ);
        end
        end
    end
    
    if east_shale==0
    else
        for j=1:2*north_shale
        for i=1:2*east_shale
            Zt(j,2*(west_shale+DX)+i,:)=Z(1,2*DX,:);
            Zbottom(j,2*(west_shale+DX)+i,1)=z2(1,2*DX,2*DZ);
        end
        end
        for j=1:2*south_shale
        for i=1:2*east_shale
            Zt(2*(north_shale+DY)+j,2*(west_shale+DX)+i,:)=Z(2*DY,2*DX,:);
            Zbottom(2*(north_shale+DY)+j,2*(west_shale+DX)+i,1)=z2(2*DY,2*DX,2*DZ);
        end
        end
    end
end

%------------------------------------------------------------------------        
    
        


Z=Zt;



%----------------------------------------------------



%start from here (last work done was on 12 FEB 2019 at 22:00 pm)




%Z1 = ones(DY3D*2,DX*2,DZ);%skip
Z3 = ones(DY3D*2,DX3D*2,DZ);% consider DX & DY instead of DX3D and DY3D
Z2 = ones(DY3D*2,DX3D*2,DZ3D);% same as above but with shale layers above and below 
%Z4 = ones(DY3D*2,DX3D*2,DZ+20);%skip


%----------------skip-----------------
%  parfor i = 1:10
%  for k = 1:9
%  Z1(i,:,k)=Z(1,:,k);
%  end
%  end
%--------------skip----------------------

 Z3(:,:,:)=Z;

%----------------skip---------------------
%  parfor j = 1:2
%  for k = 1:9
%  Z3(:,j,k)=Z1(:,1,k);
%  end
%  end
%  
%  Z3(:,1:228,:)=Z1;
%  
%  parfor j = 231:232
%  for k = 1:9
%  Z3(:,j,k)=Z1(:,228,k);
%  end
%  end
%-------------------skip---------------------- 
 
w = 100;% thickness of the shale above
 for k=1:top_shale %shale above
  Z2(:,:,k) = Z3(:,:,1)-w; %Z3(:,:,1) is the first value withen the reservoir layers
 w = w - 5; % the 10 is adjustable  
 end
 
 Z2(:,:,top_shale+1:top_shale+DZ) = Z3; % subtitute the reservoir layer into Z3 (replace 4 and 5 with the ks of the reservoir layers) 
 
 w = 5;% thickness of the each grid in the shale layer below
 for k=(top_shale+DZ+1):DZ3D %shale below
 Z2(:,:,k) = Z3(:,:,DZ)+w; %Z3(:,:,2) is the last value withen the reservoir layers
 w = w + 5;  
 end
 
 Z = Z2;
 
 
 % important to adjust the bottom most layer
 if bottom_shale==0
 else
Zbottom=Zbottom+100; % the 10 is adjustable as a value of the thickness of the shale below
 end