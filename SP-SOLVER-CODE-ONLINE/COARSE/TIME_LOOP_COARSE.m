% INPUT PARAMETERS
v = [29 30 31 32 33 34 35 36 ];% time steps
DX = 67; DY = 71; DZ = 215;% grid size


for j = 1:length(v)
    i = v(j);
formatspec = '[Uek%d,P%d,S%d,SALT%d,x,y,z,Uec%d,DX3D,DY3D,DZ3D,sigma%d,L_ec%d,L_ek%d,POROp,SATNUMp,Ute%d,TEMP%d,L_et%d] = SP_FUNCTION_COARSE(DX,DY,DZ,i);';
eval(sprintf(formatspec,i,i,i,i,i,i,i,i,i,i,i));
end
