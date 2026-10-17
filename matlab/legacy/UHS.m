function [emissivity,temperature]=UHS(ndvi,TIR,wavelength, trans, upw, dow)%, filename)

% INPUT:
% ndvi       = ndvi image
% TIR        = thermal image
% wavelength = wavelenght of TIR image in micrometers
% trans      = atmospheric transmission
% up         = upwelling radiance 
% dow        = downwelling radiance
% filename   = name of the file that temperature will be save
%
% OUTPUT:
% emissivity  = emissivity image
% temperature = temperature image

C1 = 14387.7;
C2 = 1.19104*10^8;

ndvi = double(ndvi);
TIR = double(TIR);

Nlin = size(ndvi,1);
Ncol = size(ndvi,2);

BACK = find(TIR == -2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NDVI - EMISSIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emissivity = zeros(Nlin,Ncol);

for L = 1 : Nlin
    for C = 1 : Ncol
        NDVI = ndvi(L,C);                
        %%%%%%%%%%% Emissivity from Liu and Zhang (2011)%%%%%%%%%%%%%%%%%%% 
        if (NDVI >= 0.727)
            emissivity(L,C) = 0.990;  
        end;    
        if (NDVI >= 0.157) && (NDVI <= 0.727)
            E = 1.0094 + 0.047*log(NDVI);
            emissivity(L,C) = E;  
        end;    
        if (NDVI >= -0.185) && (NDVI <= 0.157)
            emissivity(L,C) = 0.97;      
        end;    
        if (NDVI <= -0.185)
            emissivity(L,C) = 0.995;            
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end;
end;

temperature = zeros(Nlin,Ncol);
for L = 1 : Nlin
    for C = 1 : Ncol
        %%%%%%%%%%%%%%%%%%%%%%% TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A = (TIR(L,C) - upw - trans*(1 - emissivity(L,C))*dow);
        %
        B = ((wavelength^5)*A)/(trans*emissivity(L,C));
        %
        D = wavelength*log((C2/B)+1);
        %
        T = C1/D;   
        %
        temperature(L,C)= T;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end;
end;
%
temperature(BACK) = 0;
emissivity(BACK) = 0;
% %
% % imshow(temperature, []);
% % colormap(hot);
% % colorbar;
% save resultados emissivity temperature
% %
% T = temperature*100;
% %
% T = uint16(T);
% %
% filename = [filename '.tif'];
% imwrite(T,filename,'tif');
% %
end
%
%  TEMP_SCM_LOG é o nome do arquivo na area de trabalho do MATLAB
%  SCM_LOG nome do tiff
% imwrite(uint16(TEMP_SCM_LOG*100),'SCM_LOG.tiff');


