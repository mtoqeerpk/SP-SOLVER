function filecontents=readFfilezzz(filename,onlylast)
% Read the fields of an eclipse F**** of A**** format output file.
% Primarily intended for restarting eclipse using filecontents (or a
% modified version of it) to write the restartfile. Data can also be
% found from the variable filecontents.
%
% function filecontents=readFfile(filename,onlylast)
%
% filename:     name of file that is read.
% filecontents: the contents of the file is stored in this matlab
%               variable.
% onlylast:     if this variable is equal to 0, data from multiple
%               instances of the same keyword will be stored as
%               columns in filecontents.keyword.  The default
%               behavior (equivalent to onlylast==1) is to keep
%               only the last data for each keyword.
%               
%
% To produce ECLIPSE file on the required F-format, use the keyword
% 
%  RPTSCHED
%    RESTART=2 / 
% in the SCHEDULE section,
% or
%  RPTRST
%    BASIC=1 /
% in the SOLUTION section,
% and
% FMTOUT
% in the RUNSPEC section.
%
%
% Do NOT use keyword 
% UNIFOUT
% since this produces pressure and saturation data in the PRT-file.
%
% The contents is strored in the structured variable filecontens.
% filecontents.varnames and filecontens.typename contatins the list
% of variables and corresponding type that is specified in the
% eclipse file. The other entires in the filecontents variable
% contains the contents of the F file.
%
% See also writeFfile and appendFfile.
% This function is used in forweclipse.
%
% Copyright(c) RF-Rogaland Research 2002. All rights reserved.
%
% Geir Naevdal, 30.08.01.
% Speedup by AAG and GEN, Dec01-Jan02.
% AAG Feb02: Allow multiple instances of same keyword.
%            This is not guaranteed to work yet, if the number of
%            data elements changes between instances.

if nargin == 1
  onlylast=1;
end
number=1;
fid=fopen(filename,'rt');

filecontents.varnames=[];
filecontents.typename=[];

if fid<=0
  error(['Could not open file:   ',filename])
end

% Read whole file to string
s = char(fread(fid))';
%
slength=length(s);
fclose(fid);
NI=0;

cformat = ' ''%c%c%c%c%c%c%c%c''';
check = 0;
check2 = 'LGR';
check3 = 'ENDGRID';
vn = [];
% Read contents into structure

while NI < slength % there should be something left to read
    
  if strcmp(vn,check2) 
      check = 1;
      file = data;
%   elseif strcmp(vn,check3)
%       check = 0;
  elseif strcmp(vn,check2)
      check = 0;
  end
  
  NI = NI+1; % Move to first character of header line
  % Read header line
  varname = s(NI+2:NI+9);
  numentries = sscanf(s(NI+12:NI+22),'%i');
  vn = varname(~isspace(varname));
  typename = s(NI+25:NI+28);
  NI = NI+30; % Discard header line + newline character
  %
  
  if check < 1
  filecontents.varnames = [filecontents.varnames; varname];
  filecontents.typename = [filecontents.typename; typename];
  %
  if numentries > 0
    nmil = floor(numentries/1000);
    next = numentries-nmil*1000;
    switch(typename)
     case {'INTE'}
      % if typename == 'INTE'
      npos = nmil*12167+floor(next/6)*73;
      if rem(next,6) ~= 0
	npos = npos + rem(next,6)*12 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%i',numentries);
     case {'LOGI'}
      % elseif typename == 'LOGI'
      npos = nmil*3040+floor(next/25)*76;
      if rem(next,25) ~= 0
	npos = npos + rem(next,25)*3 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%s',numentries);
     case {'DOUB'}
      % elseif typename == 'DOUB'
      npos = nmil*23334+floor(next/3)*70;
      if rem(next,3) ~= 0
	npos = npos + rem(next,3)*23 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(strrep(tmp,'D','E'),'%g');
     case {'REAL'}
      % elseif typename == 'REAL'
      npos = nmil*17250+floor(next/4)*69;
      if rem(next,4) ~= 0
	npos = npos + rem(next,4)*17 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%g');
     case {'CHAR'}
      % elseif typename == 'CHAR'
      npos = nmil*11143+floor(next/7)*78;
      if rem(next,7) ~= 0
	npos = npos + rem(next,7)*11 + 1;
      end
      tmp = s(NI:NI+npos);
      data = reshape(sscanf(tmp,cformat,numentries*8),8,numentries)';
     otherwise
      % else
      error('Parse error, probably bug in readFfile.m')
    end
    NI = NI+npos;
  else
    data=[];
  end
  % The following may happen for very large or very small numbers
  % when eclipse drops the D or E for numbers with three-digit exponents.
  if size(data,1) ~= numentries && strcmp(typename,'LOGI')~=1
    error(sprintf('Did not read keyword %s correctly from file %s',...
                  varname,filename))
  end
  % If the field already exists, we try to append the data to the
  % current contents.
  vn = varname(~isspace(varname)); % The field name can not contain
                                   % spaces.
  number = number+1;
  if number >1000
      break
  end
  vn(regexp(vn,'+'))=[];
  vn(regexp(vn,'\d'))=[];
  if onlylast==0 & isfield(filecontents,vn)
    filecontents = setfield(filecontents,vn,...
		       [getfield(filecontents,vn),data]);
  else
    filecontents = setfield(filecontents,vn,data);
  end
  
  
  
  
  
  
  
  
  elseif check == 1
  formatspec = '%s.varnames = [];\n';
  eval(sprintf(formatspec,file))
  formatspec = '%s.typename = [];\n';
  eval(sprintf(formatspec,file))
  
  
  
  formatspec = '%s.varnames = [%s.varnames; varname];\n';
  eval(sprintf(formatspec,file,file))
  formatspec = '%s.typename = [%s.typename; typename];\n';
  eval(sprintf(formatspec,file,file))
  
%   filecontents.varnames = [filecontents.varnames; varname];
%   filecontents.typename = [filecontents.typename; typename];
  %
  if numentries > 0
    nmil = floor(numentries/1000);
    next = numentries-nmil*1000;
    switch(typename)
     case {'INTE'}
      % if typename == 'INTE'
      npos = nmil*12167+floor(next/6)*73;
      if rem(next,6) ~= 0
	npos = npos + rem(next,6)*12 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%i',numentries);
     case {'LOGI'}
      % elseif typename == 'LOGI'
      npos = nmil*3040+floor(next/25)*76;
      if rem(next,25) ~= 0
	npos = npos + rem(next,25)*3 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%s',numentries);
     case {'DOUB'}
      % elseif typename == 'DOUB'
      npos = nmil*23334+floor(next/3)*70;
      if rem(next,3) ~= 0
	npos = npos + rem(next,3)*23 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(strrep(tmp,'D','E'),'%g');
     case {'REAL'}
      % elseif typename == 'REAL'
      npos = nmil*17250+floor(next/4)*69;
      if rem(next,4) ~= 0
	npos = npos + rem(next,4)*17 + 1;
      end
      tmp = s(NI:NI+npos);
      data = sscanf(tmp,'%g');
     case {'CHAR'}
      % elseif typename == 'CHAR'
      npos = nmil*11143+floor(next/7)*78;
      if rem(next,7) ~= 0
	npos = npos + rem(next,7)*11 + 1;
      end
      tmp = s(NI:NI+npos);
      data = reshape(sscanf(tmp,cformat,numentries*8),8,numentries)';
     otherwise
      % else
      error('Parse error, probably bug in readFfile.m')
    end
    NI = NI+npos;
  else
    data=[];
  end
  % The following may happen for very large or very small numbers
  % when eclipse drops the D or E for numbers with three-digit exponents.
  if size(data,1) ~= numentries && strcmp(typename,'LOGI')~=1
    error(sprintf('Did not read keyword %s correctly from file %s',...
                  varname,filename))
  end
  % If the field already exists, we try to append the data to the
  % current contents.
  vn = varname(~isspace(varname)); % The field name can not contain
                                   % spaces.
  number = number+1;
  if number >1000
      break
  end
  
%   formatspec = '%s.typename = [filecontents.varnames; typename];\n';
%   eval(sprintf(formatspec,file))
  vn(regexp(vn,'+'))=[];
  vn(regexp(vn,'\d'))=[];
  if onlylast==0 & isfield(filecontents,vn)
      formatspec = '%s = setfield(%s,vn,[getfield(%s,vn),data]);\n';
      eval(sprintf(formatspec,file,file,file))       
  else
      formatspec = '%s = setfield(%s,vn,data);\n';
      eval(sprintf(formatspec,file,file))
  end 
  end
end
