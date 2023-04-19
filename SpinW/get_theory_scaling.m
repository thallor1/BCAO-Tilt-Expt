function scale_factor = get_theory_scaling(Imeas,Ierr,Itheory)    
    obs_over_err = Imeas./Ierr.^2;
    th_over_err = Itheory./Ierr.^2;
    num = sum(obs_over_err,'all');
    denom = sum(th_over_err,'all');
    scale_factor=num./denom;
    %Nan values should not contribute to this. 
end