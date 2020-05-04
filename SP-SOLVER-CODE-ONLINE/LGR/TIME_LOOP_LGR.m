%---------Input parameters--------

% Grid cells in x,y,z for the LGR
DXC = 35;
DYC = 35;
DZC = 215;

% Grid cells in x, y, z for the parent model without the added shale
DX = 67;
DY = 71;
DZ = 215;

% Grid cells in x, y, z for the parent model without the added shale
DX3D = DX;
DY3D = DY;
DZ3D = DZ+20;

% Number of LGR
num=4;

% Time steps
v=[  29 30 31 32 33 34 35 36]; 

%---------End of Input parameters--------


onlylast=1;


for j = 1:length(v)
    i = v(j);
formatspec = '[xc,yc,zc,Pc%d,Sc%d,SALTc%d,Uekc%d,Uecc%d,L_ekc%d,sigmac%d,L_ecc%d,POROc,HOSTNUM_C,TEMPc%d,Uetc%d] = SP_FUNCTION_LGR(DXC,DYC,DZC,i,onlylast,num,Uekw%d,sigma%d,P%d,L_ek%d,DX3D,DY3D,DZ3D,DX,DY,DZ,x,y,z,Uecw%d,SALT%d,L_ec%d,POROp,SATNUMp,Uetw%d,TEMP%d,L_et%d);';
eval(sprintf(formatspec,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i));
end


    