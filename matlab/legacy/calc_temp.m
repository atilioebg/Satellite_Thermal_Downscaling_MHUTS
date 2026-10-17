function [TEMPERATURA] = calc_temp (RAD,NDVI,CLASS,M,W, LD, LU, Tr)

% Downscaling AVHRR land surface temperatures 
% for improved surface urban heat
% island intensity estimation

% Using a Modified HUTS Algorithm to Downscale
% Land Surface Temperature Retrieved from Landsat-8
% Imagery: A Case Study of Xiamen City, China

% INPUTS
% RAD   : imagem radiancia termal LANDSAT8 sem correção atmosferica
% NDVI  : imagem NDVI LANDSAT8
% CLASS : imagem classificada em 3 classes (0- vegetacao, 1- area urbana, 2- agua)
% M     : equação para o calculo do Fractional Vegetation Cover (PV) 
%         (0- Gutman, 1- Carlston para veg. esp., 2- Bonafoni, 3- Chen e Ding com inv. temp.)  
% W     : water vapor content dado em g/cm^2 calculado no ENVI
% LD    : downwelling radiation calculado caso M = 3
% LU    : upwelling radiation calculada caso M = 3
% Tr    : transmittance calculado caso M = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULO DA TEMPERATURA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NDVIs = 0;

NDVIv = 0.8;

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

			if M == 0
			
				% Gutman
			
				PV = (NDVI(i,j) - NDVIs)/(NDVI(i,j) - NDVIv); 
                
            end;
			
			if M == 1
			
				% Carlston (vegetacao esparca)
				
				PV = ((NDVI(i,j) - NDVIs)/(NDVI(i,j) - NDVIv))^2; 
				
			end;
			

					if (M == 0) or (M == 1)
						
						if CLASS == 0

							emi = 0.995; % agua
                            
                        end;

						if CLASS == 1

							emi = 0.9625 + 0.00614*PV - 0.0426*(PV)^2; % super. natural
                            
                        end;
						
						if CLASS == 2

							emi = 0.9589 + 0.086*PV - 0.0671*(PV)^2; % sup. urbana
							
						end;

						Tb = K2/log((K1/RAD(i,j))+1);
						
						gama = (Tb^2)/(Br*RAD(i,j));

						delta = Tb - ((Tb^2)/Br);

						psi1 = 0.04019*(W^2) + 0.02916*W + 1.01523;

						psi2 = -0.38333*(W^2) - 1.50294*W + 0.20321;

						psi3 = 0.00918*(W^2) + 136072*W -0.27514;

						LST = gama*(((psi1*RAD(i,j) + psi2)/emi) + psi3) + delta;
						
					end;
			
			if M == 2
				
				% PV(ou fv) dado por Bonafoni, 2016 apud Stathopoulou, 2009
				
				PV = ((NDVI(i,j) - NDVI_min)/(NDVI_max - NDVI_min))^2; 
							
				emi = (1-PV)*emi_c + PV*emi_v;
				
				[LST] = TES(RAD(i,j), LAMB, emi, LD, LU, Tr);
				
			end;
			
			TEMPERATURA(i,j) = LST;

		end;
		
	end;
end		
