max_z = DZ3D;% with shale
max_x = DX3D;% with shale
max_y = DY3D;% with shale
wa=1;
wb=1;
NNCI1a=[];
NNCI2a=[];
NNCI1b=[];
NNCI2b=[];


for k=1:max_z

        for i=1:max_x
           for j=1:max_y-1
                n=(k-1)*(max_y*max_x)+(j-1)*max_x+i; %counting based on x dirction
                
                if Z(j*2,i*2,k) ~= Z(j*2+1,i*2,k)
                    if k~=max_z
                        for ww=1:max_z
                            if ww~=max_z
                                if   ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww) < Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww+1)) && (Z(j*2+1,i*2,ww+1) < Z(j*2,i*2,k+1))) ||  (Z(j*2,i*2,k) == Z(j*2+1,i*2,ww)) || (Z(j*2+1,i*2,ww+1) == Z(j*2,i*2,k+1)) || ((Z(j*2,i*2,k) > Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww+1) > Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww+1) < Z(j*2,i*2,k+1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      wa = wa+1;
                                      v = 1;
                            
                                end
                                if   ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww) < Z(j*2,i*2-1,k+1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww+1)) && (Z(j*2+1,i*2-1,ww+1) < Z(j*2,i*2-1,k+1))) ||  (Z(j*2,i*2-1,k) == Z(j*2+1,i*2-1,ww)) || (Z(j*2+1,i*2-1,ww+1) == Z(j*2,i*2-1,k+1)) || ((Z(j*2,i*2-1,k) > Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww+1) > Z(j*2,i*2-1,k+1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww+1) < Z(j*2,i*2-1,k+1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      if wa==1
                                         NNCI1a(wa) = n;
                                         NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      else
                                         if NNCI2a(wa)==NNCI2a(wa-1)
                                             NNCI1a(wa)=[];
                                             NNCI2a(wa)=[];
                                         else
                                             wa = wa+1;
                                             v = 1;
                                         end
                                       end
                                end
                            else
                                if   ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww) < Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Zbottom(j*2+1,i*2,1)) && (Zbottom(j*2+1,i*2,1) < Z(j*2,i*2,k+1))) ||  (Z(j*2,i*2,k) == Z(j*2+1,i*2,ww)) || (Zbottom(j*2+1,i*2,1) == Z(j*2,i*2,k+1)) || ((Z(j*2,i*2,k) > Z(j*2+1,i*2,ww)) && (Zbottom(j*2+1,i*2,1) > Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Zbottom(j*2+1,i*2,1) < Z(j*2,i*2,k+1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      wa = wa+1;
                                      v = 1;
                                end
                                if   ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww) < Z(j*2,i*2-1,k+1))) || ((Z(j*2,i*2-1,k) < Zbottom(j*2+1,i*2-1,1)) && (Zbottom(j*2+1,i*2-1,1) < Z(j*2,i*2-1,k+1))) ||  (Z(j*2,i*2-1,k) == Z(j*2+1,i*2-1,ww)) || (Zbottom(j*2+1,i*2-1,1) == Z(j*2,i*2-1,k+1)) || ((Z(j*2,i*2-1,k) > Z(j*2+1,i*2-1,ww)) && (Zbottom(j*2+1,i*2-1,1) > Z(j*2,i*2-1,k+1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Zbottom(j*2+1,i*2-1,1) < Z(j*2,i*2-1,k+1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                       if wa==1
                                        NNCI1a(wa) = n;
                                        NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      else
                                      if NNCI2a(wa)==NNCI2a(wa-1)
                                             NNCI1a(wa)=[];
                                             NNCI2a(wa)=[];
                                         else
                                             wa = wa+1;
                                             v = 1;
                                      end
                                      end
                                end
                                
                                
                            end
                        end
                    else
                        
                        for ww=1:max_z
                            if ww~=max_z
                                if   ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww+1)) && (Z(j*2+1,i*2,ww+1) < Zbottom(j*2,i*2,1))) || (Z(j*2,i*2,k) == Z(j*2+1,i*2,ww)) || (Z(j*2+1,i*2,ww+1) == Zbottom(j*2,i*2,1)) || ((Z(j*2,i*2,k) > Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww+1) > Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww+1) < Zbottom(j*2,i*2,1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      wa = wa+1;
                                      v = 1;
                            
                                end
                                if   ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww) < Zbottom(j*2,i*2-1,1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww+1)) && (Z(j*2+1,i*2-1,ww+1) < Zbottom(j*2,i*2-1,1))) || (Z(j*2,i*2-1,k) == Z(j*2+1,i*2-1,ww)) || (Z(j*2+1,i*2-1,ww+1) == Zbottom(j*2,i*2-1,1)) || ((Z(j*2,i*2-1,k) > Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww+1) > Zbottom(j*2,i*2-1,1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww+1) < Zbottom(j*2,i*2-1,1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      if wa==1
                                         NNCI1a(wa) = n;
                                         NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      else
                                      if NNCI2a(wa)==NNCI2a(wa-1)
                                             NNCI1a(wa)=[];
                                             NNCI2a(wa)=[];
                                         else
                                             wa = wa+1;
                                             v = 1;
                                      end
                                      end
                                end
                                
                            else
                                if   ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Z(j*2+1,i*2,ww) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Zbottom(j*2+1,i*2,1)) && (Zbottom(j*2+1,i*2,1) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) > Z(j*2+1,i*2,ww)) && (Zbottom(j*2+1,i*2,1) > Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2+1,i*2,ww)) && (Zbottom(j*2+1,i*2,1) < Zbottom(j*2,i*2,1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      wa = wa+1;
                                      v = 1;
                                end
                                if   ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Z(j*2+1,i*2-1,ww) < Zbottom(j*2,i*2-1,1))) || ((Z(j*2,i*2-1,k) < Zbottom(j*2+1,i*2-1,1)) && (Zbottom(j*2+1,i*2-1,1) < Zbottom(j*2,i*2-1,1))) || ((Z(j*2,i*2-1,k) > Z(j*2+1,i*2-1,ww)) && (Zbottom(j*2+1,i*2-1,1) > Zbottom(j*2,i*2-1,1))) || ((Z(j*2,i*2-1,k) < Z(j*2+1,i*2-1,ww)) && (Zbottom(j*2+1,i*2-1,1) < Zbottom(j*2,i*2-1,1)))
                                      NNCI1a(wa) = n;
                                      NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      if wa==1
                                          NNCI1a(wa) = n;
                                          NNCI2a(wa) = (ww-1)*(max_y*max_x)+(j)*max_x+i;
                                      else
                                      if NNCI2a(wa)==NNCI2a(wa-1)
                                             NNCI1a(wa)=[];
                                             NNCI2a(wa)=[];
                                         else
                                             wa = wa+1;
                                             v = 1;
                                      end
                                      end
                                end
                                
                            end
                        end
                    end
                end
           end
        end
end
     %------------------------------------------------------- 
     
     for k=1:max_z
        for j=1:max_y
           for i=1:max_x-1
           
                n=(k-1)*(max_y*max_x)+(j-1)*max_x+i;
                
               if Z(j*2,i*2,k) ~= Z(j*2,i*2+1,k)
                    if k~=max_z
                        for ww=1:max_z
                            if ww~=max_z
                                if   ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww) < Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww+1)) && (Z(j*2,i*2+1,ww+1) < Z(j*2,i*2,k+1))) ||  (Z(j*2,i*2,k) == Z(j*2,i*2+1,ww)) || (Z(j*2,i*2+1,ww+1) == Z(j*2,i*2,k+1)) || ((Z(j*2,i*2,k) > Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww+1) > Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww+1) < Z(j*2,i*2,k+1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      wb = wb+1;
                                      v = 1;
                            
                                end
                                if   ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww) < Z(j*2-1,i*2,k+1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww+1)) && (Z(j*2-1,i*2+1,ww+1) < Z(j*2-1,i*2,k+1))) ||  (Z(j*2-1,i*2,k) == Z(j*2-1,i*2+1,ww)) || (Z(j*2-1,i*2+1,ww+1) == Z(j*2-1,i*2,k+1)) || ((Z(j*2-1,i*2,k) > Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww+1) > Z(j*2-1,i*2,k+1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww+1) < Z(j*2-1,i*2,k+1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      if wb==1
                                          NNCI1b(wb) = n;
                                          NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      else
                                      if NNCI2b(wb)==NNCI2b(wb-1)
                                             NNCI1b(wb)=[];
                                             NNCI2b(wb)=[];
                                         else
                                             wb = wb+1;
                                             v = 1;
                                         end
                                      end
                            
                                end
                            else
                                if   ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww) < Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Zbottom(j*2,i*2+1,1)) && (Zbottom(j*2,i*2+1,1) < Z(j*2,i*2,k+1))) ||  (Z(j*2,i*2,k) == Z(j*2,i*2+1,ww)) || (Zbottom(j*2,i*2+1,1) == Z(j*2,i*2,k+1)) || ((Z(j*2,i*2,k) > Z(j*2,i*2+1,ww)) && (Zbottom(j*2,i*2+1,1) > Z(j*2,i*2,k+1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Zbottom(j*2,i*2+1,1) < Z(j*2,i*2,k+1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      wb = wb+1;
                                      v = 1;
                                end
                                if   ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww) < Z(j*2-1,i*2,k+1))) || ((Z(j*2-1,i*2,k) < Zbottom(j*2-1,i*2+1,1)) && (Zbottom(j*2-1,i*2+1,1) < Z(j*2-1,i*2,k+1))) ||  (Z(j*2-1,i*2,k) == Z(j*2-1,i*2+1,ww)) || (Zbottom(j*2-1,i*2+1,1) == Z(j*2-1,i*2,k+1)) || ((Z(j*2-1,i*2,k) > Z(j*2-1,i*2+1,ww)) && (Zbottom(j*2-1,i*2+1,1) > Z(j*2-1,i*2,k+1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Zbottom(j*2-1,i*2+1,1) < Z(j*2-1,i*2,k+1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      if wb==1
                                          NNCI1b(wb) = n;
                                          NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      else                                      
                                      if NNCI2b(wb)==NNCI2b(wb-1)
                                             NNCI1b(wb)=[];
                                             NNCI2b(wb)=[];
                                         else
                                             wb = wb+1;
                                             v = 1;
                                      end
                                      end
                                end
                            end
                        end
                    else
                        
                        for ww=1:max_z
                            if ww~=max_z
                                if   ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww+1)) && (Z(j*2,i*2+1,ww+1) < Zbottom(j*2,i*2,1))) || (Z(j*2,i*2,k) == Z(j*2,i*2+1,ww)) || (Z(j*2,i*2+1,ww+1) == Zbottom(j*2,i*2,1)) || ((Z(j*2,i*2,k) > Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww+1) > Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww+1) < Zbottom(j*2,i*2,1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      wb = wb+1;
                                      v = 1;
                            
                                end
                                if   ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww) < Zbottom(j*2-1,i*2,1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww+1)) && (Z(j*2-1,i*2+1,ww+1) < Zbottom(j*2-1,i*2,1))) || (Z(j*2-1,i*2,k) == Z(j*2-1,i*2+1,ww)) || (Z(j*2-1,i*2+1,ww+1) == Zbottom(j*2-1,i*2,1)) || ((Z(j*2-1,i*2,k) > Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww+1) > Zbottom(j*2-1,i*2,1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww+1) < Zbottom(j*2-1,i*2,1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      if wb==1
                                          NNCI1b(wb) = n;
                                          NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      else
                                      if NNCI2b(wb)==NNCI2b(wb-1)
                                             NNCI1b(wb)=[];
                                             NNCI2b(wb)=[];
                                         else
                                             wb = wb+1;
                                             v = 1;
                                      end
                                      end
                                end
                            else
                                if   ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Z(j*2,i*2+1,ww) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Zbottom(j*2,i*2+1,1)) && (Zbottom(j*2,i*2+1,1) < Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) > Z(j*2,i*2+1,ww)) && (Zbottom(j*2,i*2+1,1) > Zbottom(j*2,i*2,1))) || ((Z(j*2,i*2,k) < Z(j*2,i*2+1,ww)) && (Zbottom(j*2,i*2+1,1) < Zbottom(j*2,i*2,1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      wb = wb+1;
                                      v = 1;
                                end
                                if   ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Z(j*2-1,i*2+1,ww) < Zbottom(j*2-1,i*2,1))) || ((Z(j*2-1,i*2,k) < Zbottom(j*2-1,i*2+1,1)) && (Zbottom(j*2-1,i*2+1,1) < Zbottom(j*2-1,i*2,1))) || ((Z(j*2-1,i*2,k) > Z(j*2-1,i*2+1,ww)) && (Zbottom(j*2-1,i*2+1,1) > Zbottom(j*2-1,i*2,1))) || ((Z(j*2-1,i*2,k) < Z(j*2-1,i*2+1,ww)) && (Zbottom(j*2-1,i*2+1,1) < Zbottom(j*2-1,i*2,1)))
                                      NNCI1b(wb) = n;
                                      NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      if wb==1
                                          NNCI1b(wb) = n;
                                          NNCI2b(wb) = (ww-1)*(max_y*max_x)+(j-1)*max_x+i+1;
                                      else
                                      if NNCI2b(wb)==NNCI2b(wb-1)
                                             NNCI1b(wb)=[];
                                             NNCI2b(wb)=[];
                                         else
                                             wb = wb+1;
                                             v = 1;
                                      end
                                      end
                                end
                            end
                        end
                    end
                end         
                                   
                                
                          
    
           end
       end
     end

%-------------------------------------------------
 NNC1=[NNCI1a NNCI1b  ];
 NNC2=[ NNCI2a NNCI2b ];
 
    
     