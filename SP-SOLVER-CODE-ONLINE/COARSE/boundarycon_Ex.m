%without shale
% DX, DY, DZ
%with shale
%DX3D, DY3D, DZ3D


%--------shale------input--------
west_shale=0;
east_shale=0;
north_shale=5;
south_shale=0;
boundaryExt=zeros(DY+north_shale+south_shale,DX+west_shale+east_shale);
%------------shale-----input---------

boundaryEx = zeros(DY,DX);
% boundary(:,1) = 1;
% boundary(1,:) = 1;
% boundary(48,:) = 1;
% boundary(:,139) = 1;
[max_y, max_x]=size(boundaryEx); 

% for i=1:(max_x-1)
%      for j=1:(max_y-1)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*i,1) - z2((2*j),(2*i)+1,1);
%          
%          if un ~= 0 
%              boundaryEx(j,i+1) = 1;
%              boundaryEx(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end


% for i=1:(max_x-1)
%      for j=1:(max_y)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*i,1) - z2((2*j),(2*i)+1,1);
%          
%          if un ~= 0 
%              boundaryEx(j,i+1) = 1;
%              boundaryEx(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end

for i=1:(max_x-1)
     for j=1:(max_y)
         n=(i-1)*max_y+j;
         un = z2(2*j,2*i,1) - z2((2*j),(2*i)+1,1);
         
         if un ~= 0 
          %   boundaryEx(j,i+1) = 1;
             boundaryEx(j,i) = 1;
%              boundary(n+1) = 1;
%          else boundary(i,j) = 0;
         end
     end
end
 
%----------adjust for X and Y shale------can remove it if there is no shaqle in x and y-------------------------------

boundaryExt(north_shale+1:north_shale+DY,west_shale+1:west_shale+DX)=boundaryEx;

if north_shale==0 && south_shale==0 
    
    boundaryExt(1,west_shale+1:west_shale+DX)=boundaryEx(1,:);
    boundaryExt(DY:north_shale+DY+south_shale,west_shale+1:west_shale+DX)=boundaryEx(DY,:);

elseif north_shale==0
    boundaryExt(1,west_shale+1:west_shale+DX)=boundaryEx(1,:);

    for i=1:south_shale
       boundaryExt(north_shale+DY+i,west_shale+1:west_shale+DX)=boundaryEx(DY,:);
    end


elseif south_shale==0

    for i=1:north_shale
       boundaryExt(i,west_shale+1:west_shale+DX)=boundaryEx(1,:);
    end
    boundaryExt(north_shale+DY:north_shale+DY+south_shale,west_shale+1:west_shale+DX)=boundaryEx(DY,:);

else

    for i=1:north_shale
       boundaryExt(i,west_shale+1:west_shale+DX)=boundaryEx(1,:);
    end
    for i=1:south_shale
       boundaryExt(north_shale+DY+i,west_shale+1:west_shale+DX)=boundaryEx(DY,:);
    end

end 

boundaryEx=boundaryExt;









