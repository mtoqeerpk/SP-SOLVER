function [Uek,P,S,SALT,x,y,z,Uec,DX3D,DY3D,DZ3D,sigma,L_ec,L_ek,PORO3D,SATNUM_SHALE,Ute,TEMP,L_et] = SP_FUNCTION_COARSE( DX, DY, DZ, ntime)

% -------------------------READ GRID FILES------------------------------------
    EGRID ='FILE NAME.FEGRID'; %FILE NAME FROM ECLIPSE
    
    %[F]=readFfile(EGRID);%readFfile to be considered for files with NO LGR
    [F]=readFfilezz(EGRID);%readFfilezz to be considered when reading compiled file of parent/LGR
    
    % x,y,z coordinates   
    
    COORD = F.COORD;
    COORD1= reshape(COORD,6,[]);
    XCORD = COORD1(1,1:DX+1);
    ylength = size(COORD1);
    YCORD = COORD1(2,1:DX+1:ylength(2));  
    ZCORN = F.ZCORN;
    z1 = reshape(ZCORN,2*DX,2*DY,2*DZ);
    parfor ip = 1:2*DZ
        z2(:,:,ip) = z1(:,:,ip)';
     end;
    z3 = z2(:,[1:2:(2*DX-1),2*DX],:);
    z4 = z3([1:2:(2*DY-1),2*DY],:,:);
    z5 = z4(:,:,[1:2:(2*DZ-1),2*DZ]);
    ZCORD = z5(1,1,1:DZ+1);

  
    XCORD = double(XCORD);
    YCORD = double(YCORD);
    ZCORD = double(ZCORD);
    
   % Converting ft to m 
    x = XCORD*0.3048;
    y = YCORD*0.3048;
    z = ZCORD*0.3048;
  
    
    ATN = F.ACTNUM; 
    ACTNUM1 = reshape(ATN, DX, DY, DZ); 

    parfor ip = 1:DZ
        ACTNUM(:,:,ip) = ACTNUM1(:,:,ip)';
    end
    
% -------------------------END OF READ GRID FILES------------------------------------ 

% -------------------------READ F FILES------------------------------------

  % Read F Files at n (crtain Time step)

    if ntime <= 9
    FFILE =['FILE NAME.F000', int2str(ntime)];%FILE NAME FROM ECLIPSE
    %G = readFfile(FFILE);%readFfile to be considered for files with NO LGR
    G = readFfilezzz(FFILE);%readFfilezzz to be considered when reading compiled file of parent/LGR
    
elseif (9<ntime) && (ntime<=99)
    FFILE =['FILE NAME.F00', int2str(ntime)];%FILE NAME FROM ECLIPSE
    %G = readFfile(FFILE);%readFfile to be considered for files with NO LGR
    G = readFfilezzz(FFILE);%readFfilezzz to be considered when reading compiled file of parent/LGR
    
elseif (99<ntime) && (ntime<=999)
    FFILE =['FILE NAME.F0', int2str(ntime)];%FILE NAME FROM ECLIPSE
    %G = readFfile(FFILE);%readFfile to be considered for files with NO LGR
    G = readFfilezzz(FFILE);%readFfilezzz to be considered when reading compiled file of parent/LGR
    
elseif (999<ntime) && (ntime<=9999)
    FFILE =['FILE NAME.F', int2str(ntime)];%FILE NAME FROM ECLIPSE
    %G = readFfile(FFILE);%readFfile to be considered for files with NO LGR
    G = readFfilezzz(FFILE);%readFfilezzz to be considered when reading compiled file of parent/LGR
    
    end
  

% Read SW, Pressure, and Salt concentration and TEMPERATURE into 3D Matrix.   
SWAT = G.SWAT(:);
SALT = G.SALT(:); 
PRES = G.WAT_POTN(:); 
TEMP = G.TEMP(:); 

% SW Matrix
jp = 1;
    for r = 1:length(ATN)        
        if ATN(r) == 0            
            saturation(r) = 0;            
        else
            saturation(r) = SWAT(jp);
            jp = jp+1;        
        end        
    end
    
X = reshape(saturation, DX, DY, DZ);


    parfor ip = 1:DZ
        
        SWAT3D(:,:,ip) = X(:,:,ip)';
    end 

    
clear j
jp=1;

% Pressure Matrix
    for r = 1:length(ATN)
        if ATN(r) == 0
        pressure(r) = 0;
        else
            pressure(r) = PRES(jp);
        jp = jp+1;
        end  
    end
    
X = reshape(pressure, DX, DY, DZ);

    parfor ip = 1:DZ
        PRES3D(:,:,ip) = X(:,:,ip)';
    end
    
clear j
jp=1;

% Salt Matrix
     for r = 1:length(ATN)
         if ATN(r) == 0
         salt(r) = 40.97; %max Salt in the reservoir
         else
             salt(r) = SALT(jp);
         jp = jp+1;
         end  
     end
     
 X = reshape(salt, DX, DY, DZ);
 
     for ip = 1:DZ
         SALT3D(:,:,ip) = X(:,:,ip)';
     end

 
clear j
jp=1;

% Temperature Matrix
     for r = 1:length(ATN)
         if ATN(r) == 0
         temperature(r) = 343; % max temperature in the reservoir
         else
             temperature(r) = TEMP(jp);
         jp = jp+1;
         end  
     end
     
 X = reshape(temperature, DX, DY, DZ);
 
     for ip = 1:DZ
         TEMP3D(:,:,ip) = X(:,:,ip)';
     end
     
  clear j
% -------------------------END OF READ F FILES------------------------------------     

% ----------------------ADD POROSITY AND SATURATION REGION AND FINAL DEFINITION OF VARIABLES---------------------------
% Porosity reading

load ('poro_2');%POROSITY LOADED FROM ECLIPSE
PORO1=poro_2;

% porosity Matrix
    
X = reshape(PORO1, DX, DY, DZ);

    parfor ip = 1:DZ
        PORO3D(:,:,ip) = X(:,:,ip)';
    end
    
clear j

% SATNUM TO IDENTIFY THE SATURATION REGIONS
load('SATNUM_1.mat')% SATNUM TO BE LOADED FROM ECLIPSE

is=1;
for k=1:DZ
    for j=1:DY
        for i=1:DX
            SATNUM(j,i,k)=SATNUM_1(is);
            is=is+1;
        end
    end
end

% Final definition of the main variables and unit conversion.

PRES3D = double(PRES3D);
SWAT3D = double(SWAT3D);
SALT3D = double(SALT3D); 
TEMP3D = double(TEMP3D); 

 S = SWAT3D;
 P = PRES3D*6894;      % convert from psi to Pa
 SALT = SALT3D.*2.852801./58.44; % convert from lbm/stb to mole/litre
 TEMP =(TEMP3D-32)*5/9+273.15; % convert from F to K
 
% ---------------------END OF ADD POROSITY AND SATURATION REGION AND FINAL DEFINITION OF VARIABLES----------------------

%---------Identify Boundaries for each Cell-----------------------------------

% %need to enter manually the NORTH,SOUTH,EAST,WEST shales. not required if
% %there is no faults.
% 
%  boundarycon_Ex; 
%  boundarycon_Wx; 
%  boundarycon_Ey; 
%  boundarycon_Wy;  
 
%---------End of Identify Boundaries for each Cell-----------------------------------

% ----------------------Adding Shale Layers--------------------------------

%-----------------------------CASE A---------------------------------------
 % CASE A: shale layers above and below only (Z-direction). 
 
    DX3D = DX;
    DY3D = DY;
    DZ3D = DZ+20;  % 20 cells extra. can be modified to any other number.
    
 % creat 3D matrix for S,P,SALT,POROSITY and temperature to account for the expanded
 % dimension 
 
    S3D = ones(DY3D,DX3D,DZ3D);  % Fully saturate with water.
    P3D = 6894*2100*ones(DY3D,DX3D,DZ3D);  % initial pressure.
    POROS3D = 0.3*ones(DY3D,DX3D,DZ3D);  % assumed porosity value for shale
    SALT3D = 40.97*(2.852801/58.44)*ones(DY3D,DX3D,DZ3D);  % formation water salinity
    SATNUM_SHALE = zeros(DY3D,DX3D,DZ3D);
    TEMP3D = 343*ones(DY3D,DX3D,DZ3D); % reservoir temerature 
    
 % Assign values of S, P, SALT, POROSITY and temperature for the non shale layers.
 % Can be modified based on total number of shale layers considered.
    S3D(:,:,11:DZ3D-10)=S; 
    P3D(:,:,11:DZ3D-10)=P;
    SALT3D(:,:,11:DZ3D-10)=SALT;
    POROS3D(:,:,11:DZ3D-10)=PORO3D;
    SATNUM_SHALE(:,:,11:DZ3D-10)= SATNUM;
    TEMP3D(:,:,11:DZ3D-10)=TEMP;
    
 % Define the X,Y,Z coordinate matrix including the shale layers. 
    y3D = ones(1,DY3D+1);
    z3D = ones(1,DZ3D+1);
    x3D = ones(1,DX3D+1);

    %if no shale in x and y
     y3D(:)=y; % Simple replacment
     x3D(:)=x; % Simple replacment
       
 % 70 m of shales above and below. Any thickness can be considered with
 % considering adjusting the thickness of each shale layer (below)
 
    z3D(1:11) = [(z(1)-70) (z(1)-63) (z(1)-56) (z(1)-49) (z(1)-42) (z(1)-35) (z(1)-28) (z(1)-21) (z(1)-14) (z(1)-7) z(1)];           
    z3D(11:DZ3D-9) = z;
    z3D(DZ+12:DZ3D+1) = [(z(end)+7) (z(end)+14) (z(end)+21) (z(end)+28) (z(end)+35) (z(end)+42) (z(end)+49) (z(end)+56) (z(end)+63) (z(end)+70)];
%-----------------------------END OF CASE A---------------------------------------

%-----------------------------CASE B---------------------------------------
% % CASE B: shale layers in x and y and z .    
% 
% %      80 shale layers in z (40 A&B each 10m) & 60 shale layes in x
% %      (30 E&W each 3.3528m) & 6 shale layers in y (3 E&W each 33.528 m)
% %      Can be modified based on any other number of shale layers considered.
% 
%     DX3D = DX+60;
%     DY3D = DY+6;
%     DZ3D = DZ+80;  
% 
% % creat 3D matrix for S,P,SALT,POROSITY and temperature to account for the expanded
%  % dimension 
%  
%     S3D = ones(DY3D,DX3D,DZ3D);  % Fully saturate with water.
%     P3D = 6894*2100*ones(DY3D,DX3D,DZ3D);  % initial pressure
%     POROS3D = 0.3*ones(DY3D,DX3D,DZ3D);  % assumed porosity value for shale
%     SALT3D = 40.97*(2.852801/58.44)*ones(DY3D,DX3D,DZ3D);  % formation water salinity
%     SATNUM_SHALE = zeros(DY3D,DX3D,DZ3D);
%     TEMP3D = 343*ones(DY3D,DX3D,DZ3D); % reservoir temerature 

%     S3D(4:DY3D-3,31:DX3D-30,41:DZ3D-40)=S;
%     P3D(4:DY3D-3,31:DX3D-30,41:DZ3D-40)=P;
%     SALT3D(4:DY3D-3,31:DX3D-30,41:DZ3D-40)=SALT;
%     POROS3D(4:DY3D-3,31:DX3D-30,41:DZ3D-40)=PORO3D; 
%     SATNUM_SHALE(4:DY3D-3,31:DX3D-30,41:DZ3D-40)= SATNUM;
%     TEMP3D(4:DY3D-3,31:DX3D-30,41:DZ3D-40)=TEMP;

% % Define the X,Y,Z coordinate matrix including the shale layers. 
%     y3D = ones(1,DY3D+1);
%     z3D = ones(1,DZ3D+1);
%     x3D = ones(1,DX3D+1);
%     
% 
% %     %80 shale layers in z (40 A&B each 10m) & 60 shale layes in x
% %     %(30 E&W each 3.3528m) & 6 shale layers in y (3 E&W each 33.528 m)
% %     % Can be modified based any number of shale layers considered
% 
%      z3D(1:41) = [(z(1)-400) (z(1)-390) (z(1)-380) (z(1)-370) (z(1)-360) (z(1)-350) (z(1)-340) (z(1)-330) (z(1)-320) (z(1)-310) (z(1)-300) (z(1)-290) (z(1)-280) (z(1)-270) (z(1)-260) (z(1)-250) (z(1)-240) (z(1)-230) (z(1)-220) (z(1)-210) (z(1)-200) (z(1)-190) (z(1)-180) (z(1)-170) (z(1)-160) (z(1)-150) (z(1)-140) (z(1)-130) (z(1)-120) (z(1)-110) (z(1)-100) (z(1)-90) (z(1)-80) (z(1)-70) (z(1)-60) (z(1)-50) (z(1)-40) (z(1)-30) (z(1)-20) (z(1)-10) z(1)];           
%      z3D(41:DZ3D-39) = z;
%      z3D(DZ+42:DZ3D+1) = [(z(end)+10) (z(end)+20) (z(end)+30) (z(end)+40) (z(end)+50) (z(end)+60) (z(end)+70) (z(end)+80) (z(end)+90) (z(end)+100) (z(end)+110) (z(end)+120) (z(end)+130) (z(end)+140) (z(end)+150) (z(end)+160) (z(end)+170) (z(end)+180) (z(end)+190) (z(end)+200) (z(end)+210) (z(end)+220) (z(end)+230) (z(end)+240) (z(end)+250) (z(end)+260) (z(end)+270) (z(end)+280) (z(end)+290) (z(end)+300) (z(end)+310) (z(end)+320) (z(end)+330) (z(end)+340) (z(end)+350) (z(end)+360) (z(end)+370) (z(end)+380) (z(end)+390) (z(end)+400)];
%  
%      x=x+100.584;
%      x3D(1:31)= [(x(1)-90.584) (x(1)-87.584) (x(1)-84.584) (x(1)-81.584) (x(1)-78.584) (x(1)-75.584) (x(1)-72.584) (x(1)-69.584) (x(1)-66.584) (x(1)-63.584) (x(1)-60.584) (x(1)-57.584) (x(1)-54.584) (x(1)-51.584) (x(1)-48.584) (x(1)-45.584) (x(1)-42.584) (x(1)-39.584) (x(1)-36.584) (x(1)-33.584) (x(1)-30.584) (x(1)-27.584) (x(1)-24.584) (x(1)-21.584) (x(1)-18.584) (x(1)-15.584) (x(1)-12.584) (x(1)-9.584) (x(1)-6.584) (x(1)-3.3528) x(1)];
%      x3D(31:DX3D-29)= x;
%      x3D(DX+32:DX3D+1)= [(x(end)+3.3528) (x(end)+6.584) (x(end)+9.584) (x(end)+12.584) (x(end)+15.584) (x(end)+18.584) (x(end)+21.584) (x(end)+24.584) (x(end)+27.584) (x(end)+30.584) (x(end)+33.584) (x(end)+36.584) (x(end)+39.584) (x(end)+42.584) (x(end)+45.584) (x(end)+48.584) (x(end)+51.584) (x(end)+54.584) (x(end)+57.584) (x(end)+60.584) (x(end)+63.584) (x(end)+66.584) (x(end)+69.584) (x(end)+72.584) (x(end)+75.584) (x(end)+78.584) (x(end)+81.584) (x(end)+84.584) (x(end)+87.584) (x(end)+90.584)];  
%      
%      y=y+100.584;
%      y3D(1:4)= [(y(1)-100.584) (y(1)-67.056) (y(1)-33.528) y(1)];
%      y3D(4:DY3D-2)= y;
%      y3D(DY+5:DY3D+1)=[(y(end)+33.528) (y(end)+67.056) (y(end)+100.584)];

%-----------------------------END OF CASE B---------------------------------------
 
 % Final definition of the variables including the shale layers.
    P = P3D;
    S = S3D;
    SALT = SALT3D;
    PORO3D = POROS3D;
    TEMP = TEMP3D;
    x = x3D;
    y = y3D;
    z = z3D;
    
 % ----------------------End of Adding Shale Layers--------------------------------  
 
% clean the WORKSPACE
    clearvars -except x y z S P SALT TEMP F G DX DY DZ DX3D DY3D DZ3D boundary NNC NNCI1a NNCI2a NNCI1b NNCI2b PORO3D ACTNUM ntime boundaryEx boundaryEy boundaryWx boundaryWy z2 SATNUM

%-------------------------ELECTRICAL PROPERTIES----------------------------------
% ELECTRICAL PROPERTIES

% Define Swirr, Sor
    
Swirr = 0.2*ones(DY3D,DX3D,DZ3D); % user defined
Sor =0.2*ones(DY3D,DX3D,DZ3D); % user defined

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
    
    tna = reshape(tna, DY3D, DX3D, DZ3D);
    
    %Exclusion end member at Swirr  
    C_ee = (-0.0861.*TEMP).*10^-3;   
    %Diffusion end member at (1-Sor)
    C_ed = (-0.0861.*TEMP.*(2*tna-1)).*10^-3; 
    m_ec = 3;
    Cr_ec = Cec_relative(Swn, m_ec);
    
    % C_ec ecxlusion-diffusion affect
    EPS=0.4; % Shale exclusion efficiency- user defined
    C_ec(:,:,:) = EPS*(C_ee(:,:,:)-C_ed(:,:,:))+C_ed(:,:,:);
    
    % CASE A. Shale above and below only. modified based on layers
    % considered for Case A
    C_ec(:,:,11:DZ3D-10) = (1-Cr_ec(:,:,11:DZ3D-10)).*(C_ed(:,:,11:DZ3D-10)-C_ee(:,:,11:DZ3D-10))+C_ee(:,:,11:DZ3D-10);
    
%     % CASE B. Shale in x, y and z. modified based on layers
%     % considered for Case B
%     C_ec(4:DY3D-3,31:DX3D-30,41:DZ3D-40) = (1-Cr_ec(4:DY3D-3,31:DX3D-30,41:DZ3D-40)).*(C_ed(4:DY3D-3,31:DX3D-30,41:DZ3D-40)-C_ee(4:DY3D-3,31:DX3D-30,41:DZ3D-40))+C_ee(4:DY3D-3,31:DX3D-30,41:DZ3D-40);



% C_te calculation. unit is V/(K).
    %Exclusion end member at Swirr  
    C_eet = ((-1.984/10)*(log10(SALT))+(5.953/10)).*10^-3;    
    %Diffusion end member at (1-Sor)
    C_edt = ((-1.984/10)*(2.*tna-1).*(log10(SALT))+1.059.*tna-(5.673/10)).*10^-3; 
    m_te = 3;
    Cr_te = Cte_relative(Swn, m_te);
    
    % C_ec ecxlusion-diffusion affect
    EPS=0.4; % Shale exclusion efficiency- user defined
    C_te(:,:,:) = EPS*(C_eet(:,:,:)-C_edt(:,:,:))+C_edt(:,:,:);
    
    % CASE A. Shale above and below only. modified based on layers
    % considered for Case A
    C_te(:,:,11:DZ3D-10) = (1-Cr_te(:,:,11:DZ3D-10)).*(C_edt(:,:,11:DZ3D-10)-C_eet(:,:,11:DZ3D-10))+C_eet(:,:,11:DZ3D-10);
    
%     % CASE B. Shale in x, y and z. modified based on layers
%     % considered for Case B
%     C_te(4:DY3D-3,31:DX3D-30,41:DZ3D-40) = (1-Cr_te(4:DY3D-3,31:DX3D-30,41:DZ3D-40)).*(C_edt(4:DY3D-3,31:DX3D-30,41:DZ3D-40)-C_eet(4:DY3D-3,31:DX3D-30,41:DZ3D-40))+C_eet(4:DY3D-3,31:DX3D-30,41:DZ3D-40);


% Water conductivity Calculation. 
  
    for r = 1:length(SALT(:))              
      sigma1(r) = (5.6 + 0.27*(TEMP(r)-273.15)- 1.5*(10^-4)*((TEMP(r)-273.15)^2))*SALT(r)-((2.36+0.099*(TEMP(r)-273.15))/(1+0.214*((SALT(r))^1/2)))*((SALT(r))^(3/2));          
    end

  sigma1 = reshape( sigma1, DY3D, DX3D, DZ3D);
 
% Saturated rock conductivity.

 sigma = sigma_Swa(sigma1, S, PORO3D);

%Calculating the coupling terms
    L_ec = sigma.*C_ec;
    L_ek = sigma.*C_ek.*Cr_ek;
    L_te = sigma.*C_te;
      
% Set L_ek = 0 for shale layers assuming no flow across shale layers.
 
ATN4D = zeros(DY3D,DX3D,DZ3D);
        
%     % CASE B. Shale in x, y and z. modified based on layers
%     % considered for Case A
      ATN4D(:,:,11:10+DZ) = ACTNUM;
      
%     % CASE B. Shale in x, y and z. modified based on layers
%     % considered for Case B   
%     ATN4D(4:3+DY,31:30+DX,41:40+DZ) = ATN3D;


 % set Lek =0    
    for k=1:DZ3D
        for j=1:DX3D
            for i=1:DY3D
                n=(k-1)*(DY3D*DX3D)+(j-1)*DY3D+i;
                
           
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

 SP_SIGMA_3D;
%SP_FAULT_SIGMA_3D;
 
 SP_Lek_3D_JIA_forcedzero;
%SP_FAULT_Lek_3D_JIA_forcedzero;
    
 SP_Lec_3D_JIA_HIGH_CF;
%SP_FAULT_Lec_3D_JIA_HIGH_CF;

 SP_Lte_3D_JIA;
%SP_FAULT_Lte_3D_JIA;

%Set any erroneous values in sparse matrices to zero
    A(isnan(A)) = 0;
    B(isnan(B)) = 0;    
    C(isnan(C)) = 0;      
    D(isnan(D)) = 0; 

%---Generate NNC List------------    
% % Extend the Z2 to shale layer to be able to generate NNC list
 
%  Extenting_the_coordinates_shale;% need to adjust manually for the surrounding shales   
%  NNC_shale;% no thing to chang 
%  NNC_LIST_FAULT;%no thing to change
%---End of Generate NNC List------------ 


%------NNC adjusting---------
%  SP_NNC1_Ex_SIGMA;
%  SP_NNC1_Ey_SIGMA;
%  SP_NNC2_Wx_SIGMA;
%  SP_NNC2_Wy_SIGMA;

%  SP_NNC1_Ex_LEK_JIA_FORCEDZERO;
%  SP_NNC1_Ey_LEK_JIA_FORCEDZERO;
%  SP_NNC2_Wx_LEK_JIA_FORCEDZERO; 
%  SP_NNC2_Wy_LEK_JIA_FORCEDZERO;

%  SP_NNC1_Ex_LEC_JIA_HIGH_CF; 
%  SP_NNC1_Ey_LEC_JIA_HIGH_CF; 
%  SP_NNC2_Wx_LEC_JIA_HIGH_CF; 
%  SP_NNC2_Wy_LEC_JIA_HIGH_CF; 

%  SP_NNC1_Ex_LTE_JIA; 
%  SP_NNC1_Ey_LTE_JIA; 
%  SP_NNC2_Wx_LTE_JIA; 
%  SP_NNC2_Wy_LTE_JIA; 

% %Set any erroneous values in sparse matrices to zero
%     A(isnan(A)) = 0;
%     B(isnan(B)) = 0;    
%     C(isnan(C)) = 0;  
%     D(isnan(D)) = 0; 
%------End of NNC adjusting---------

    
    
%Implement direct solver to calculate x from Ax = B*P
    %[Uek,flag,relres,iter,resvec] = bicgstab(-A,B*P(:),1.0e-9,10000,[]);
    Uek = mldivide(-A,B*P(:));  % folwing lines for masseges-later
    con=sprintf(' Uek model solved by brute force for TS%d',ntime);
    disp(con)
    
%Implement direct solver to calculate x from Ax = C*Cf    
   %[Uec,flag,relres,iter,resvec] = bicgstab(-A,C*log(SALT(:)),1.0e-12,10000,[]);
   Uec = mldivide(-A,C*log(SALT(:)));  % folwing lines for masseges-later
    con=sprintf(' Uec model solved by brute force for TS%d',ntime);
    disp(con)     
    
%Implement direct solver to calculate x from Ax = D*T    
   %[Ute,flag,relres,iter,resvec] = bicgstab(-A,D*TEMP(:),1.0e-12,10000,[]);
   Ute = mldivide(-A,D*TEMP(:));  % folwing lines for masseges-later
    con=sprintf(' Ute model solved by brute force for TS%d',ntime);
    disp(con)
    
%Output maximum errors from direct solver
    errek = max(abs(-A*Uek-B*P(:)));
    disp(errek)
    
    errec = max(abs(-A*Uec-C*log(SALT(:))));
    disp(errec) 
    
    errec = max(abs(-A*Ute-C*TEMP(:)));
    disp(errec) 
    
Uek=reshape(Uek,DY3D,DX3D,DZ3D);
Uec=reshape(Uec,DY3D,DX3D,DZ3D);
Ute=reshape(Ute,DY3D,DX3D,DZ3D);
    
end 


    