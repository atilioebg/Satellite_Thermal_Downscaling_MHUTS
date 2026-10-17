function [EMI,PV] = sentinel_pve (NDVI)

% Calcula e emissividade e o PV a partir do NDVI da imagem SENTINEL

NDVI = double(NDVI);
NDVI_min = min(min(NDVI));
NDVI_max = max(max(NDVI));

    for i = 1 : size(NDVI,1)
        for j = 1 : size(NDVI,2)	

            ndvi = NDVI(i,j);                

            % PV(ou fv) dado por Bonafoni, 2016 apud Stathopoulou, 2009				
            PV(i,j) = ((NDVI(i,j) - NDVI_min)/(NDVI_max - NDVI_min))^2; 

            % Emissivity from Liu and Zhang (2011) 
            if (ndvi >= 0.727)                
                emi = 0.990;                
            end;    

            if (ndvi >= 0.157) && (ndvi <= 0.727)                
                E = 1.0094 + 0.047*log(ndvi);                    
                emi = E;                  
            end;    

            if (ndvi >= -0.185) && (ndvi <= 0.157)                
                emi = 0.97;               
            end;    

            if (ndvi <= -0.185)                
                emi = 0.995;                 
            end;
            EMI(i,j) = emi;

        end;
    end; 
end

