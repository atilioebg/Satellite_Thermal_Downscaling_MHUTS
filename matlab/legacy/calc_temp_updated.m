function [TEMPERATURA] = calc_temp_updated (RAD, NDVI, CLASS, W, LD, LU, Tr, MT, ME)

% Downscaling AVHRR land surface temperatures for improved surface urban 
% heat island intensity estimation
%
% Using a Modified HUTS Algorithm to Downscale Land Surface Temperature 
% Retrieved from Landsat-8 Imagery: A Case Study of Xiamen City, China
%
% INPUTS
% RAD   : imagem radiancia termal LANDSAT8 sem correção atmosferica
% NDVI  : imagem NDVI LANDSAT8
% CLASS : imagem classificada em 3 classes (1- solo, 2- vegetacao, 3 - area urbana, 4- agua)
% W     : water vapor content dado em g/cm^2 calculado no ENVI (em cm)
% LD    : downwelling radiation calculado caso MT = 1
% LU    : upwelling radiation calculada caso MT = 1
% Tr    : transmittance calculado caso MT = 1
% MT    : como calcula temperatua (0- SCM, 1- RCM)
% ME    : como calcula emissividade (00- LOG-NDVI, 11- PV-NDVI)

% O SC (MT=0) com a emissividade a partir do LOG (ME=00) pareceu o mais plausível
RAD = double(RAD);

NDVI = double(NDVI);

K1 = 774.89; % W / m^2 sr micro

K2 = 1321.08; % K

Br = 1324; % K 

NDVI_min = min(min(NDVI));

NDVI_max = max(max(NDVI));

emi_c = 0.99;

emi_v = 0.94; 

% LANDSAT 8 central wav. http://www.gisagmaps.com/landsat-8-sentinel-2-bands/
LAMB =  10.895; % wavelenght in micrometer (without cientific notation(x10^y))

	for i = 1 : size(RAD,1)

		for j = 1 : size(RAD,2)
            
            %%%%%%%%%%%%%%%%%%% CALCULO DA EMISSIVIDADE %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if ME == 00
                
                ndvi = NDVI(i,j);                
                
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
                
            end;
            
            if ME == 11
                            
                % PV(ou fv) dado por Bonafoni, 2016 apud Stathopoulou, 2009
				
				PV = ((NDVI(i,j) - NDVI_min)/(NDVI_max - NDVI_min))^2; 
                
                classe = CLASS(i,j);
            
				if classe == 4
                    
                    % agua
                    emi = 0.995;
                
                else
                    
                    % demais alvos
                    emi = (1-PV)*emi_c + PV*emi_v;
                    
                end;
                
            end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %
            %   
            %%%%%%%%%%%%%%%%%%% CALCULO DA TEMPERATURA %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if MT == 0

                % Single Channel Method (SCM)
                Tb = K2/log((K1/RAD(i,j))+1);

                gama = (Tb^2)/(Br*RAD(i,j));

                delta = Tb - ((Tb^2)/Br);

                psi1 = 0.04019*(W^2) + 0.02916*W + 1.01523;

                psi2 = -0.38333*(W^2) - 1.50294*W + 0.20321;

                psi3 = 0.00918*(W^2) + 1.36072*W -0.27514;

                LST = gama*(((psi1*RAD(i,j) + psi2)/emi) + psi3) + delta;

            end;

            if MT == 1

                % Reference Channel Method (RCM)
                [LST] = TES(RAD(i,j), LAMB, emi, LD, LU, Tr);

            end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			TEMPERATURA(i,j) = LST;

		end;
		
	end;
    
end		

