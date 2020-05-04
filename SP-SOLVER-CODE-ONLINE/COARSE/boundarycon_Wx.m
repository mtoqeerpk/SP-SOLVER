%without shale
% DX, DY, DZ
%with shale
%DX3D, DY3D, DZ3D

%--------shale------input--------
west_shale=0;
east_shale=0;
north_shale=5;
south_shale=0;
boundaryWxt=zeros(DY+north_shale+south_shale,DX+west_shale+east_shale);
%------------shale-----input---------

boundaryWx = zeros(DY,DX);
% boundary(:,1) = 1;
% boundary(1,:) = 1;
% boundary(48,:) = 1;
% boundary(:,139) = 1;
[max_y, max_x]=size(boundaryWx); 


% for i=2:(max_x)
%      for j=2:(max_y)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*(i-1),1) - z2((2*j),(2*i)-1,1);
%          
%          if un ~= 0 
%              boundaryWx(j,i-1) = 1;
%              boundaryWx(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end

% for i=2:(max_x)
%      for j=1:(max_y)
%          n=(i-1)*max_y+j;
%          un = z2(2*j,2*(i-1),1) - z2((2*j),(2*i)-1,1);
%          
%          if un ~= 0 
%              boundaryWx(j,i-1) = 1;
%              boundaryWx(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end


for i=2:(max_x)
     for j=1:(max_y)
         n=(i-1)*max_y+j;
         un = z2(2*j,2*(i-1),1) - z2((2*j),(2*i)-1,1);
         
         if un ~= 0 
             %boundaryWx(j,i-1) = 1;
             boundaryWx(j,i) = 1;
%              boundary(n+1) = 1;
%          else boundary(i,j) = 0;
         end
     end
end
 

%----------adjust for X and Y shale------can remove it if there is no shaqle in x and y-------------------------------

boundaryWxt(north_shale+1:north_shale+DY,west_shale+1:west_shale+DX)=boundaryWx;

if north_shale==0 && south_shale==0 
    
    boundaryWxt(1,west_shale+1:west_shale+DX)=boundaryWx(1,:);
    boundaryWxt(DY:north_shale+DY+south_shale,west_shale+1:west_shale+DX)=boundaryWx(DY,:);

elseif north_shale==0
    boundaryWxt(1,west_shale+1:west_shale+DX)=boundaryWx(1,:);

    for i=1:south_shale
       boundaryWxt(north_shale+DY+i,west_shale+1:west_shale+DX)=boundaryWx(DY,:);
    end


elseif south_shale==0

    for i=1:north_shale
       boundaryWxt(i,west_shale+1:west_shale+DX)=boundaryWx(1,:);
    end
    boundaryWxt(north_shale+DY:north_shale+DY+south_shale,west_shale+1:west_shale+DX)=boundaryWx(DY,:);

else

    for i=1:north_shale
       boundaryWxt(i,west_shale+1:west_shale+DX)=boundaryWx(1,:);
    end
    for i=1:south_shale
       boundaryWxt(north_shale+DY+i,west_shale+1:west_shale+DX)=boundaryWx(DY,:);
    end

end 

boundaryWx=boundaryWxt;

