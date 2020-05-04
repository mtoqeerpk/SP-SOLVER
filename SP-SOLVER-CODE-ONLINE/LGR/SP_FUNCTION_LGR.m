function [x,y,z,P,S,SALT,Uekc,Uecc,L_ek,sigma,L_ec,PORO3D,HOSTNUM_C,TEMP,Uetc]= SP_FUNCTION_LGR(DXC,DYC,DZC,ntime,onlylast, num, Uekn,sigman,Pn,L_ekn,DX3D,DY3D,DZ3D,DX,DY,DZ,xp,yp,zp,Uecn,SALTn,L_ecn,POROp,SATNUMp,Uetn,TEMPn,Letn)

% -------------------------READ GRID FILES------------------------------------
filename ='FILE NAME.FEGRID'; %FILE NAME FROM ECLIPSE
readLGRGRID;

if num <= 9

	formatspec = 'COORD= LGR%d.COORD;\n';
	eval(sprintf(formatspec,num));
	formatspec = 'ZCORN = LGR%d.ZCORN;\n';
	eval(sprintf(formatspec,num));

    formatspec = 'HOSTNUM = LGR%d.HOSTNUM;\n';
	eval(sprintf(formatspec,num));
   
elseif (9<num) && (num<=999)

	formatspec = 'COORD= LGR%d.COORD;\n';
	eval(sprintf(formatspec,num));
	formatspec = 'ZCORN = LGR%d.ZCORN;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'HOSTNUM = LGR%d.HOSTNUM;\n';
	eval(sprintf(formatspec,num));
end

filename ='FILE NAME.FEGRID'; %FILE NAME FROM ECLIPSE
readLGRGRID;

if num <= 9
	formatspec = 'ATN = LGR%d.ACTNUM;\n';
	eval(sprintf(formatspec,num));  
elseif (9<num) && (num<=999)
	formatspec = 'ATN = LGR%d.ACTNUM;\n';
	eval(sprintf(formatspec,num));
end

ACTNUM1 = reshape(ATN, DXC, DYC, DZC);

    parfor ip = 1:DZ
        
        ACTNUM(:,:,ip) = ACTNUM1(:,:,ip)';
    end

% x,y,z coordinates

    COORD1= reshape(COORD,6,[]);
    XCORD = COORD1(1,1:DXC+1);
    ylength = size(COORD1);
    YCORD = COORD1(2,1:DXC+1:ylength(2));  
    z1 = reshape(ZCORN,2*DXC,2*DYC,2*DZC);
    for ip = 1:2*DZC
        z2(:,:,ip) = z1(:,:,ip)';
    end
    z3 = z2(:,[1:2:(2*DXC-1),2*DXC],:);
    z4 = z3([1:2:(2*DYC-1),2*DYC],:,:);
    z5 = z4(:,:,[1:2:(2*DZC-1),2*DZC]);
    ZCORD = z5(1,1,1:DZC+1);
    
    XCORD = double(XCORD);
    YCORD = double(YCORD);
    ZCORD = double(ZCORD);
    
    % Converting ft to m 
    x = XCORD*0.3048;
    y = YCORD*0.3048;
    z = ZCORD*0.3048;
    

% define the host cells for each child cell.
    jp = 1;
    for r = 1:length(ATN)        
         HN(r) = HOSTNUM(jp);            
            jp = jp+1; 
    end
    
X = reshape(HN, DXC, DYC, DZC);


    parfor ip = 1:DZC
        
        HN3D(:,:,ip) = X(:,:,ip)';
    end  
clear j
  
% -------------------------END OF READ GRID FILES------------------------------------ 

% -------------------------READ F FILES------------------------------------

  % Read F Files at n (crtain Time step)
  
if ntime <= 9
    filename =['%FILE NAME FROM ECLIPSE.F000', int2str(ntime)]; %FILE NAME FROM ECLIPSE
    readLGRFFILE;
    
elseif (9<ntime) && (ntime<=99)
    filename =['%FILE NAME FROM ECLIPSE.F00', int2str(ntime)]; %FILE NAME FROM ECLIPSE
    readLGRFFILE;
    
elseif (99<ntime) && (ntime<=999)
    filename =['%FILE NAME FROM ECLIPSE.F0', int2str(ntime)]; %FILE NAME FROM ECLIPSE
    readLGRFFILE;

elseif (999<ntime) && (ntime<=9999)
    filename =['%FILE NAME FROM ECLIPSE.F', int2str(ntime)]; %FILE NAME FROM ECLIPSE
    readLGRFFILE;
end
  
% Read SW, Pressure, and Salt concentration and TEMPERATURE into 3D Matrix. 

 if num <= 9
	formatspec = 'SWAT = LGR%d.SWAT;\n';
	eval(sprintf(formatspec,num));
	formatspec = 'PRES= LGR%d.WAT_POTN;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'HOSTNUM= LGR%d.HOSTNUM;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'SALT= LGR%d.SALT;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'TEMP= LGR%d.TEMP;\n';
	eval(sprintf(formatspec,num));
    
elseif (9<num) && (num<=99)
	formatspec = 'SWAT = LGR%d.SWAT;\n';
	eval(sprintf(formatspec,num));
	formatspec = 'PRES= LGR%d.WAT_POTN(:);\n';
	eval(sprintf(formatspec,num));
    formatspec = 'HOSTNUM= LGR%d.HOSTNUM;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'SALT= LGR%d.SALT;\n';
	eval(sprintf(formatspec,num));
    formatspec = 'TEMP= LGR%d.TEMP;\n';
	eval(sprintf(formatspec,num));
 end


% SW Matrix
jp = 1;
    for r = 1:length(ATN)        
        if ATN(r) == 0            
            saturation(r) = 0;            
        else saturation(r) = SWAT(jp);            
            jp = jp+1;        
        end        
    end
    
X = reshape(saturation, DXC, DYC, DZC);


    parfor ip = 1:DZC
        
        SWAT3D(:,:,ip) = X(:,:,ip)';
    end  
clear j

jp=1;

% Pressure Matrix
    for r = 1:length(ATN)
        if ATN(r) == 0
        pressure(r) = 0;
        else pressure(r) = PRES(jp);
        jp = jp+1;
        end  
    end
    
X = reshape(pressure, DXC, DYC, DZC);

    parfor ip = 1:DZC
        PRES3D(:,:,ip) = X(:,:,ip)';
    end
    
clear j

jp=1;

% Salt Matrix
     for r = 1:length(ATN)
         if ATN(r) == 0
         salt(r) = 40.97; %max Salt in the reservoir
         else salt(r) = SALT(jp);
         jp = jp+1;
         end  
     end
     
 X = reshape(salt, DXC, DYC, DZC);
 
     for ip = 1:DZC
         SALT3D(:,:,ip) = X(:,:,ip)';
     end
     
clear j

jp=1;

% Temperature Matrix
     for r = 1:length(ATN)
         if ATN(r) == 0
         temperature(r) = 343; %max temperature in the reservoir
         else temperature(r) = TEMP(jp);
         jp = jp+1;
         end  
     end
     
 X = reshape(temperature, DXC, DYC, DZC);
 
     for ip = 1:DZC
         TEMP3D(:,:,ip) = X(:,:,ip)';
     end
 clear j    


% Final definition of the main variables and unit conversion. 

PRES3D = double(PRES3D);
SWAT3D = double(SWAT3D);
SALT3D = double(SALT3D);
TEMP3D = double(TEMP3D);

 S = SWAT3D;
 P = PRES3D*6894;      % convert from psi to Pa
 SALT = SALT3D.*2.852801./58.44; % convert from lbm/stb to mole/litre
 HOSTNUMi = HN3D;
 TEMP =(TEMP3D-32)*5/9+273.15; % convert from F to K
 
% -------------------------END OF READ F FILES------------------------------------ 

DX3DC = DXC;
DY3DC = DYC;
DZ3DC = DZC;

% ----------------------ADD POROSITY AND SATURATION REGION---------------------------

HOST_CELLS_FINAL; 

% porosity 
PORO3D = zeros(DY3DC,DX3DC,DZ3DC);
for k=1:DZ3DC
    for i=1:DX3DC
        for j=1:DY3DC
            PORO3D(j,i,k)=POROp(HOSTNUM_C(j,i,k));
        end
    end
end


% SATNUM TO IDENTIFY THE SATURATION REGIONS
SATNUMc=zeros(DY3DC,DX3DC,DZ3DC);
for k=1:DZ3DC
    for i=1:DX3DC
        for j=1:DY3DC
            SATNUMc(j,i,k)=SATNUMp(HOSTNUM_C(j,i,k));
        end
    end
end

% ----------------------END OF ADD POROSITY AND SATURATION REGION---------------------------

%-------------------------ELECTRICAL PROPERTIES----------------------------------
% ELECTRICAL PROPERTIES

% Define Swirr, Sor

Swirr = 0.2*ones(DY3DC,DX3DC,DZ3DC); % user defined
Sor =0.2*ones(DY3DC,DX3DC,DZ3DC); % user defined

 
 % Normalized water saturation
 Swn = (S-Swirr)./(1.0-Swirr-Sor);
 Swn(Swn(:)>1)=1;   
 Swn(Swn(:)<0)=0;
 
 
% C_ek calculation. unit is V/Pa.

C_ek = (-1.36*(SALT.^-0.9123)).*10^-9;

m_ek = 0.6;
Cr_ek = Cek_relative(Swn, m_ek);
 

% C_ec calculation. unit is V/(M).

% tna calculation as a matrix.
for r = 1:length(SALT(:))        
        if SALT(r) < 0.09            
            tna(r) = 0.39;
   
       elseif SALT(r) > 0.09
           tna(r) = 0.366-(0.0212*(log10(SALT(r)))); 
           
        end        
    end
    
    tna = reshape(tna, DY3DC, DX3DC, DZ3DC);
    
    %Exclusion end member at Swirr  
    C_ee = (-0.0861.*TEMP).*10^-3;   
    %Diffusion end member at (1-Sor)
    C_ed = (-0.0861.*TEMP.*(2*tna-1)).*10^-3;  
    m_ec = 3;
    Cr_ec = Cec_relative(Swn, m_ec);
 
    %at intermediat saturation
    C_ec(:,:,:) = (1-Cr_ec(:,:,:)).*(C_ed(:,:,:)-C_ee(:,:,:))+C_ee(:,:,:); 
    

% C_te calculation. unit is V/(K).

    %Exclusion end member at Swirr  
    C_eet = ((-1.984/10)*(log10(SALT))+(5.953/10)).*10^-3;    
    %Diffusion end member at (1-Sor)
    C_edt = ((-1.984/10)*(2.*tna-1).*(log10(SALT))+1.059.*tna-(5.673/10)).*10^-3; 
    m_te = 3;
    Cr_te = Cte_relative(Swn, m_te);

    %at intermediat saturation
    C_te(:,:,:) = (1-Cr_te(:,:,:)).*(C_edt(:,:,:)-C_eet(:,:,:))+C_eet(:,:,:); 

 % Water conductivity Calculation.
 
 for r = 1:length(SALT(:))        
   sigma1(r) = (5.6 + 0.27*(TEMP(r)-273.15)- 1.5*(10^-4)*((TEMP(r)-273.15)^2))*SALT(r)-((2.36+0.099*(TEMP(r)-273.15))/(1+0.214*((SALT(r))^1/2)))*((SALT(r))^(3/2));                
 end
     
  sigma1 = reshape( sigma1, DY3DC, DX3DC, DZ3DC);
  
 
 % Saturated rock conductivity.

 sigma = sigma_Swa(sigma1, S, PORO3D);

 
 %Calculating the coupling terms
    L_ec = sigma.*C_ec;
    L_ek = sigma.*C_ek.*Cr_ek;
    L_te = sigma.*C_te;
    
 
% Set L_ek = 0  non active cells.
 
ATN4D = zeros(DY3DC,DX3DC,DZ3DC);         
ATN4D(:,:,:) = ACTNUM;
 
% set Lek =0    
    for k=1:DZ3DC
        for j=1:DX3DC
            for i=1:DY3DC
                n=(k-1)*(DY3DC*DX3DC)+(j-1)*DY3DC+i;
                
           
                if ATN4D(n) == 0
                        L_ek(n) = 0;
                        
                else
                    L_ek(n) = L_ek(n);
                    
                end          
                
            end
        end
    end 
    
%-------------------------END OF ELECTRICAL PROPERTIES----------------------------------

%-------------------------SP SOLUTION----------------------------------

% Run the SIGMA_3d, L_ek, and L_ec loops for creating to solution matrix.
 
 SP_SIGMA_3D_LGR;

 SP_Lek_3D_JIA_forcedzero_LGR;
 SP_SIGMA_3D_STUek;
 SP_SIGMA_3D_STUp_forcedzero;

 SP_Lec_3D_JIA_HIGH_CF_LGR;
 SP_SIGMA_3D_STUec;
 SP_SIGMA_3D_STUsalt_HIGH_CF;
 
 SP_Lte_3D_JIA_LGR;
 SP_SIGMA_3D_STUte;
 SP_SIGMA_3D_STUtemp;

 
 %Set any erroneous values in sparse matrices to zero
    A(isnan(A)) = 0;
    B(isnan(B)) = 0;
    C(isnan(C)) = 0;
    STUek(isnan(STUek)) = 0;
    STUp(isnan(STUp)) = 0;
    STUec(isnan(STUec)) = 0;
    STUsalt(isnan(STUsalt)) = 0;
    STUte(isnan(STUte)) = 0;
    STUtemp(isnan(STUtemp)) = 0;
 

 %Implement direct solver to calculate x from Ax = B*P + STUp - STUek
    %[Uekc,flag,relres,iter,resvec] = bicgstab(A,B*P(:)+STUp(:)-STUek(:),1.0e-6,10000,[]);
    Uekc = mldivide(A,B*P(:)+STUp(:)-STUek(:));  % folwing lines for masseges-later
    con=sprintf('Uek model solved by brute force for TS%d',ntime);
    disp(con)  
    
%Implement direct solver to calculate x from Ax = C*Cf + STUsalt - STUec   
    %[Uecc,flag,relres,iter,resvec] = bicgstab(A,C*SALT(:)+STUsalt(:)-STUec(:),1.0e-9,10000,[]);
    Uecc = mldivide(A,C*log(SALT(:))+STUsalt(:)-STUec(:)); 
    con=sprintf('Uec model solved by brute force for TS%d',ntime);
    disp(con)  

 %Implement direct solver to calculate x from Ax = C*TEMP + STUtemp - STUte
    %[Uekc,flag,relres,iter,resvec] = bicgstab(A,D*TEMP(:)+STUtemp(:)-STUte(:),1.0e-6,10000,[]);
    Utec = mldivide(A,D*TEMP(:)+STUtemp(:)-STUte(:));  % folwing lines for masseges-later
    con=sprintf('Ute model solved by brute force for TS%d',ntime);
    disp(con) 
    
   
%Output maximum errors from direct solver
    errek = max(abs(A*Uekc-B*P(:)-STUp(:)+STUek(:)));
    disp(errek)
    
    errek = max(abs(A*Uecc-C*log(SALT(:))-STUsalt(:)+STUec(:)));
    disp(errek)
    
    errek = max(abs(A*Utec-D*TEMP(:)-STUtemp(:)+STUte(:)));
    disp(errek)
    
    
    %Reshape solution matrices into 3D format including the shale layers
    Uekc = reshape(Uekc,DY3DC,DX3DC,DZ3DC);
    Uecc = reshape(Uecc,DY3DC,DX3DC,DZ3DC);
    Utec = reshape(Utec,DY3DC,DX3DC,DZ3DC);
    
end
