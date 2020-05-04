% %parent model
 %DX3D=4;
 %DY3D=3;
 %DZ3D=4;
 %DX=4;
 %DY=3;
 %DZ=2;

% %child model
% DX3DC=3;
% DY3DC=3;
% DZ3DC=9;

FL = DX3D*DY3D;

% HSX = (DX3D-DX)/2; % half of the shale layers
% HSY = (DY3D-DY)/2; % half of the shale layers
% HSZ = (DZ3D-DZ)/2; % half of the shale layers

HSX = 0; % half of the shale layers
HSY = 0; % half of the shale layers
HSZ = 0; % half of the shale layers



for n = 1:length(NNC1)
      
    if NNC1(n) > FL
    
        k = 1 + floor(NNC1(n)/FL)+ HSZ;
        if (floor(NNC1(n)/FL)-(NNC1(n)/FL))==0
            k=k-1;
        end
    
        j = 1 + floor(rem(NNC1(n),(FL*(k-1-HSZ)))/DX3D)+HSY;
        if (rem(NNC1(n),(FL*(k-1-HSZ))))==0
            j=DY3D+HSY;
        else
            if (rem(NNC1(n),(FL*(k-1-HSZ)))/DX3D) - floor(rem(NNC1(n),(FL*(k-1-HSZ)))/DX3D) == 0
            j = j - 1 ;
            end
        end
        
        i = rem(NNC1(n),(FL*(k-1-HSZ)))-(DX3D*(j-1-HSY))+HSX;
        if i <= 0
            i = DX3D + HSX;
        end
    
    else
    
         k = 1 + HSZ;
   
        j = 1 + floor(NNC1(n)/DX3D) + HSY;
        if (NNC1(n)/DX3D) - floor((NNC1(n)/DX3D)) == 0
            j = j - 1 ;
        end
   
        i = NNC1(n) - (DX3D*(j-1-HSY))+ HSX; % this is the coorect one 
        %i = HOSTNUMi(n) - (DX*(j-1))+ HSX; 
%         if j == 0
%             j = DX + HSX;
%         end
    
    end
    
    NNCC1(n,:) = [ j i k ];

end

for n = 1:length(NNC2)
      
    if NNC2(n) > FL
    
        k = 1 + floor(NNC2(n)/FL)+ HSZ;
        if (floor(NNC2(n)/FL)-(NNC2(n)/FL))==0
            k=k-1;
        end
    
        j = 1 + floor(rem(NNC2(n),(FL*(k-1-HSZ)))/DX3D)+HSY;
        if (rem(NNC2(n),(FL*(k-1-HSZ))))==0
            j=DY3D+HSY;
        else
            if (rem(NNC2(n),(FL*(k-1-HSZ)))/DX3D) - floor(rem(NNC2(n),(FL*(k-1-HSZ)))/DX3D) == 0
            j = j - 1 ;
            end
        end
        
        i = rem(NNC2(n),(FL*(k-1-HSZ)))-(DX3D*(j-1-HSY))+HSX;
        if i <= 0
            i = DX3D + HSX;
        end
    
    else
    
         k = 1 + HSZ;
   
        j = 1 + floor(NNC2(n)/DX3D) + HSY;
        if (NNC2(n)/DX3D) - floor((NNC2(n)/DX3D)) == 0
            j = j - 1 ;
        end
   
        i = NNC2(n) - (DX3D*(j-1-HSY))+ HSX; % this is the coorect one 
        %i = HOSTNUMi(n) - (DX*(j-1))+ HSX; 
%         if j == 0
%             j = DX + HSX;
%         end
    
    end
    
    NNCC2(n,:) = [ j i k ];

end