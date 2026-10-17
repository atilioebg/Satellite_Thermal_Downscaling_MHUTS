function [ALBEDO_IMG] = albedo_sentinel (REFLEC,M,I)

% Capability of Sentinel-2 data for estimating maximum 
% evapotranspiration and irrigation requirements for 
% tomato crop in Central Italy 

% Cross-Comparison of Albedo Products for Glacier
% Surfaces Derived from Airborne and Satellite
% (Sentinel-2 and Landsat 8) Optical Data

% INPUTS
% REFLEC: imagem SENTINEL de reflectancia com atm corrigida 
% M     : equação para o cáculo do albedo (1- Knap, 2- LIANG, 3- D'Urso)

% OUTPUT
% ALBEDO_IMG: imagem albedo
% NDVI      : imagem NDVI
% NDBI      : imagem NDBI
% NDWI      : imagem NDWI
 
%%%%%%%%%%%%%%%%%%%%%% CALCULO DO ALBEDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REFLEC = double(REFLEC)/10000;

for i = 1 : size(REFLEC,1)

	for j = 1 : size(REFLEC,2)
	
	B2 = REFLEC(i,j,1);
	
	B3 = REFLEC(i,j,2);
	
	B4 = REFLEC(i,j,3);
	
	B5 = REFLEC(i,j,4);
	
	B6 = REFLEC(i,j,5);
	
	B7 = REFLEC(i,j,6);
	
	B8 = REFLEC(i,j,7);
	
	B11 = REFLEC(i,j,8);
	
	B12 = REFLEC(i,j,9);
	
	if M == 1
		
		% KNAP
		
		albedo = 0.726*B3 - 0.0322*(B3)^2 - 0.015*B8 + 0.581*(B8)^2;
    
    end;
		
	if M == 2
		
		% LIANG

		albedo = 0.356*B2 + 0.130*B4 + 0.373*B8 + 0.085*B11 + 0.072*B12 - 0.0018;
        
    end;
    
	if M == 3

		% D'Urso and Calera, 2006 <============ ESTE PARECE O MAIS CORRETO
		
		albedo = 0.1324*B2 + 0.1269*B3 + 0.1051*B4 + 0.0971*B5 + 0.089*B6 + 0.0818*B7 + 0.0722*B8 + 0.0167*B11 + 0.0002*B12;
		
	end;
	
	ALBEDO_IMG(i,j) = albedo;


% %%%%%%%%%%%%%%%%%%% CALCULO DOS INDICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 	if I == 4
% 
% 		NDVI(i,j) = (B8 - B4)/(B8 + B4);
% 		
% 		NDBI(i,j) = (B11 - B8)/(B11 + B8);
% 			
% 		NDWI(i,j) = (B3- B8)/(B3 + B8);
% 			
% 	end;
%     
%     end;
    
end;

end


