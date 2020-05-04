%without shale
% DX, DY, DZ
%with shale
%DX3D, DY3D, DZ3D

%--------shale------input--------
west_shale=0;
east_shale=0;
north_shale=5;
south_shale=0;
boundaryEyt=zeros(DY+north_shale+south_shale,DX+west_shale+east_shale);

%------------shale-----input---------

boundaryEy = zeros(DY,DX);
% boundary(:,1) = 1;
% boundary(1,:) = 1;
% boundary(48,:) = 1;
% boundary(:,139) = 1;
[max_y, max_x]=size(boundaryEy);

% for i=1:(max_x-1)
%      for j=1:(max_y-1)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*i,1) - z2((2*j)+1,(2*i),1);
%          
%          if un ~= 0 
%              boundaryEy(j+1,i) = 1;
%              boundaryEy(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
% end

% for i=1:(max_x)
%      for j=1:(max_y-1)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*i,1) - z2((2*j)+1,(2*i),1);
%          
%          if un ~= 0 
%              boundaryEy(j+1,i) = 1;
%              boundaryEy(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
% end

for i=1:(max_x)
     for j=1:(max_y-1)
         n=(i-1)*max_y+j;
         un = z2(2*j,2*i,1) - z2((2*j)+1,(2*i),1);
         
         if un ~= 0 
             %boundaryEy(j+1,i) = 1;
             boundaryEy(j,i) = 1;
%              boundary(n+1) = 1;
%          else boundary(i,j) = 0;
         end
     end
end
 
 
%----------adjust for X and Y shale------can remove it if there is no shaqle in x and y-------------------------------

boundaryEyt(north_shale+1:north_shale+DY,west_shale+1:west_shale+DX)=boundaryEy;

if west_shale==0 && east_shale==0 
    
    boundaryEyt(north_shale+1:north_shale+DY,1)=boundaryEy(:,1);
    boundaryEyt(north_shale+1:north_shale+DY,DX:west_shale+DX+east_shale)=boundaryEy(:,DX);

elseif west_shale==0
    boundaryEyt(north_shale+1:north_shale+DY,1)=boundaryEy(:,1);

    for i=1:east_shale
       boundaryEyt(north_shale+1:north_shale+DY,west_shale+DX+i)=boundaryEy(:,DX);
    end


elseif east_shale==0

    for i=1:west_shale
       boundaryEyt(north_shale+1:north_shale+DY,i)=boundaryEy(:,1);
    end
    boundaryEyt(north_shale+1:north_shale+DY,west_shale+DX:west_shale+DX+east_shale)=boundaryEy(:,DX);

else

    for i=1:west_shale
       boundaryEyt(north_shale+1:north_shale+DY,i)=boundaryEy(:,1);
    end
    for i=1:east_shale
       boundaryEyt(north_shale+1:north_shale+DY,west_shale+DX+i)=boundaryEy(:,DX);
    end

end 

boundaryEy=boundaryEyt;

 
 