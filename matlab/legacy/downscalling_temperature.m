function [TRES, ESTD] = downscalling_temperature(SREFL,LRAD, SNDVI, LNDVI, WPC)


[SALB, ~, ~, ~] = albedo_sentinel (SREFL,3,4);

[SEMI,SPV] = sentinel_pve (SNDVI);

[TEMPERATURA] = calc_temp_updated (LRAD, LNDVI, 000, WPC, 000, 000, 000, 0, 00);

[TRES, ESTD] = MMQ4((TEMPERATURA*100), SEMI, SALB, SPV);


end


