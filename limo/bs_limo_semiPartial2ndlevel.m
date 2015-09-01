function [] = bs_limo_semiPartial2ndlevel(LIMO)
cd([LIMO.dir '/groupLimo'])
load('LIMO.mat')
for k = 1:5
    data = bs_limoGatherData('paths',LIMO.data.data_dir,'whatToLoad','semi_partial_coef','predictor',k);
    
    
    LIMO.design.bootstrap = 1000;
    bs_limo_random_robust(1,data,200+k,LIMO.design.bootstrap,LIMO.design.tfce);
end

