% %parent model
% DX3D=3;
% DY3D=3;
% DZ3D=3;
% DX=3;
% DY=3;
% DZ=3;

% %child model
% DX3DC=3;
% DY3DC=3;
% DZ3DC=9;


FL = DX*DY;

HSX = (DX3D-DX)/2; % half of the shale layers
HSY = (DY3D-DY)/2; % half of the shale layers
HSZ = (DZ3D-DZ)/2; % half of the shale layers


for n = 1:length(HOSTNUMi(:))
    
    if HOSTNUMi(n) > FL
    
        k = 1 + floor(HOSTNUMi(n)/FL)+ HSZ;
        if (floor(HOSTNUMi(n)/FL)-(HOSTNUMi(n)/FL))==0
            k=k-1;
        end
    
        j = 1 + floor(rem(HOSTNUMi(n),(FL*(k-1-HSZ)))/DX)+HSY;
        if (rem(HOSTNUMi(n),(FL*(k-1-HSZ))))==0
            j=DY+HSY;
        else
            if (rem(HOSTNUMi(n),(FL*(k-1-HSZ)))/DX) - floor(rem(HOSTNUMi(n),(FL*(k-1-HSZ)))/DX) == 0
            j = j - 1 ;
            end
        end
        
        i = rem(HOSTNUMi(n),(FL*(k-1-HSZ)))-(DX*(j-1-HSY))+HSX;
        if i <= 0
            i = DX + HSX;
        end
    
    else
    
        k = 1 + HSZ;
   
        j = 1 + floor(HOSTNUMi(n)/DX) + HSY;
        if (HOSTNUMi(n)/DX) - floor((HOSTNUMi(n)/DX)) == 0
            j = j - 1 ;
        end
   
        i = HOSTNUMi(n) - (DX*(j-1-HSY))+ HSX; % this is the coorect one 
        %i = HOSTNUMi(n) - (DX*(j-1))+ HSX; 
%         if j == 0
%             j = DX + HSX;
%         end
    
    end
    
    nn=(k-1)*(DY3D*DX3D)+(i-1)*DY3D+j;
    HOSTNUM_C(n) = nn;% corrected for Counting (C)
    
    %------------------------------------------------------
    %added to fix the LGR study
    iparent(n)=i;
    jparent(n)=j;
    kparent(n)=k;
    %-----------------------------------------------------
end
 

HOSTNUM_C = reshape (HOSTNUM_C,DY3DC,DX3DC,DZ3DC);
    %------------------------------------------------------
    %added to fix the LGR study
    iparent = reshape (iparent,DY3DC,DX3DC,DZ3DC);
    jparent = reshape (jparent,DY3DC,DX3DC,DZ3DC);
    kparent = reshape (kparent,DY3DC,DX3DC,DZ3DC);
    %-----------------------------------------------------

% HOSTNUM_C = reshape (HOSTNUM_C,DYC,DXC,DZC);
%     %------------------------------------------------------
%     %added to fix the LGR study
%     iparent = reshape (iparent,DYC,DXC,DZC);
%     jparent = reshape (jparent,DYC,DXC,DZC);
%     kparent = reshape (kparent,DYC,DXC,DZC);
%     %-----------------------------------------------------
%     HOSTNUM_C3D = ones(DY3DC,DX3DC,DZ3DC);
%     HOSTNUM_C3D(:,:,21:DZ3DC-20)=HOSTNUM_C;
%     
%     iparent3D = ones(DY3DC,DX3DC,DZ3DC);
%     iparent3D(:,:,21:DZ3DC-20)=iparent;
%     
%     jparent3D = ones(DY3DC,DX3DC,DZ3DC);
%     jparent3D(:,:,21:DZ3DC-20)=jparent;
%     
%     kparent3D = ones(DY3DC,DX3DC,DZ3DC);
%     kparent3D(:,:,21:DZ3DC-20)=kparent;
%     
%     for k=1:20
%         kk=21-k;
%         HOSTNUM_C3D(:,:,kk)=HOSTNUM_C(:,:,1)-k*DY3D*DX3D;
%         iparent3D(:,:,kk)=iparent(:,:,1);
%         jparent3D(:,:,kk)=jparent(:,:,1);
%         kparent3D(:,:,kk)=kparent(:,:,1)-k;
%     end
%     
%     for k=57:76
%         kk=k-56;
%         HOSTNUM_C3D(:,:,k)=HOSTNUM_C(:,:,36)+kk*DY3D*DX3D;
%         iparent3D(:,:,k)=iparent(:,:,36);
%         jparent3D(:,:,k)=jparent(:,:,36);
%         kparent3D(:,:,k)=kparent(:,:,36)+kk;
%     end
%     
%     HOSTNUM_C = HOSTNUM_C3D;
%     iparent = iparent3D;
%     jparent = jparent3D;
%     kparent = kparent3D;
%         
%         
%         
%         
%     
%     
% 
