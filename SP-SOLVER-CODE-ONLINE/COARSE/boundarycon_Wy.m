%without shale
% DX, DY, DZ
%with shale
%DX3D, DY3D, DZ3D

%--------shale------input--------
west_shale=0;
east_shale=0;
north_shale=5;
south_shale=0;
boundaryWyt=zeros(DY+north_shale+south_shale,DX+west_shale+east_shale);
%------------shale-----input---------

boundaryWy = zeros(DY,DX);
% boundary(:,1) = 1;
% boundary(1,:) = 1;
% boundary(48,:) = 1;
% boundary(:,139) = 1;
[max_y, max_x]=size(boundaryWy);

% for i=2:(max_x)
%      for j=2:(max_y)
%          n=(i-1)*max_y+j;
%          un = z2(2*(j-1),2*i,1) - z2((2*j)-1,(2*i),1);
%          
%          if un ~= 0 
%              boundaryWy(j-1,i) = 1;
%              boundaryWy(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end

% for i=1:(max_x)
%      for j=2:(max_y)
%          n=(i-1)*max_y+j;
%          un = z2(2*(j-1),2*i,1) - z2((2*j)-1,(2*i),1);
%          
%          if un ~= 0 
%              boundaryWy(j-1,i) = 1;
%              boundaryWy(j,i) = 1;
% %              boundary(n+1) = 1;
% %          else boundary(i,j) = 0;
%          end
%      end
%  end

for i=1:(max_x)
     for j=2:(max_y)
         n=(i-1)*max_y+j;
         un = z2(2*(j-1),2*i,1) - z2((2*j)-1,(2*i),1);
         
         if un ~= 0 
             %boundaryWy(j-1,i) = 1;
             boundaryWy(j,i) = 1;
%              boundary(n+1) = 1;
%          else boundary(i,j) = 0;
         end
     end
end
 
%----------adjust for X and Y shale------can remove it if there is no shaqle in x and y-------------------------------

boundaryWyt(north_shale+1:north_shale+DY,west_shale+1:west_shale+DX)=boundaryWy;

if west_shale==0 && east_shale==0 
    
    boundaryWyt(north_shale+1:north_shale+DY,1)=boundaryWy(:,1);
    boundaryWyt(north_shale+1:north_shale+DY,DX:west_shale+DX+east_shale)=boundaryWy(:,DX);

elseif west_shale==0
    boundaryWyt(north_shale+1:north_shale+DY,1)=boundaryWy(:,1);

    for i=1:east_shale
       boundaryWyt(north_shale+1:north_shale+DY,west_shale+DX+i)=boundaryWy(:,DX);
    end


elseif east_shale==0

    for i=1:west_shale
       boundaryWyt(north_shale+1:north_shale+DY,i)=boundaryWy(:,1);
    end
    boundaryWyt(north_shale+1:north_shale+DY,west_shale+DX:west_shale+DX+east_shale)=boundaryWy(:,DX);

else

    for i=1:west_shale
       boundaryWyt(north_shale+1:north_shale+DY,i)=boundaryWy(:,1);
    end
    for i=1:east_shale
       boundaryWyt(north_shale+1:north_shale+DY,west_shale+DX+i)=boundaryWy(:,DX);
    end

end 

boundaryWy=boundaryWyt;
