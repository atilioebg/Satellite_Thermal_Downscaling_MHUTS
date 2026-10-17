function [T,ESTD] = MMQ4(LST, EMI, ALB, PV)


% li =1;
% ci =1;
% 
% lf = size(EMI,1); 
% cf = size(EMI,2);



li =640; %lst 160
ci =500; %lst 125

lf =960; %lst 240
cf =760; %lst 190

% LST = double(LST)/100;

Y = uint16(zeros(size(EMI,1),size(EMI,2)));

I = 1;
	for i = 2 : 3 : size(EMI,1)-1
        J = 1;
		for j = 2 : 3: size(EMI,2)-1
					
			
			Y(i-1,j-1)= LST(I,J);
			Y(i-1,j  )= LST(I,J);
			Y(i-1,j+1)= LST(I,J);
			Y(i  ,j-1)= LST(I,J);
			Y(i  ,j  )= LST(I,J);
			Y(i  ,j+1)= LST(I,J);
			Y(i+1,j-1)= LST(I,J);
			Y(i+1,j  )= LST(I,J);
			Y(i+1,j+1)=	LST(I,J);
					
			J = J + 1;
		
		end;
		I = I + 1;
	end;
    
    % correção da borda vertical direita
    I = 1;
    fim = size(EMI,2);
    for i = 1 : 3 : size(EMI,1)-2
        Y([i:i+2],[(fim-1) fim]) = LST(I,size(LST,2));
        I = I + 1;
	end;
    
    % correção da borda horizontal inferior
    J = 0;
    fim = size(EMI,1);
    for i = 1 : 3 : size(EMI,2)-2
        J = J + 1;
        Y([(fim-1) fim],[i:i+2]) = LST(size(LST,1),J);
	end;
    
    % correção do canto inferior direito
    FL = size(EMI,1);
    FC = size(EMI,2);
    Y([FL-1 FL],[FC-1 FC]) = LST(size(LST,1),size(LST,2));
    
    EMI = EMI([li:lf],[ci:cf]);
    ALB = ALB([li:lf],[ci:cf]);
    PV  =  PV([li:lf],[ci:cf]);
    Y   =   Y([li:lf],[ci:cf]);
    
    % reformatação da imagem
    a = reshape(ALB,size(EMI,1)*size(EMI,2),1);
    e = reshape(EMI,size(EMI,1)*size(EMI,2),1);
    p = reshape(PV, size(EMI,1)*size(EMI,2),1);
    y = double(reshape(Y, size(EMI,1)*size(EMI,2),1)/100);
    
    % ajuste por MMQ
    FIT = fitlm([p e a] ,y,'poly333');
    % previsão usando o MMQ
    t = predict(FIT,[p e a]);
    
    T = reshape(t, size(EMI,1), size(EMI,2));
    
    
    MEDIA  = mean(t); 
    DESVIO = std(t); 
    RMSE   = sqrt(mean((y-t).^2));
    MAX    = max(t);
    MIN    = min(t);

    ESTD = [MEDIA DESVIO RMSE MAX MIN]';
		
end