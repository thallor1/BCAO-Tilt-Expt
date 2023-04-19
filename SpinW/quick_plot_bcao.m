function quick_plot_bcao(J1x,J1z,J3x,J3z,gab)
    % Simple script to nicely plot bcao spec assuming the main script exists already. 

    %First import data
    %close('all')
    gz=2.0;
    f_ei6_3T = 'Data files/H00 vs E B=3T_slice.dat';
    f_ei6_2T = 'Data files/H00 vs E B=2T_slice.dat';
    f_ei6_4T = 'Data files/H00 vs E B=4T_slice.dat';
    f_ei6_5T = 'Data files/H00 vs E B=5T_slice.dat';
    f_ei27_3T = 'Data files/Ei27_3T_ascii_h00.dat';
    %Store the relevant information in an object.
    BCAO_6_2T = BCAO_meas_spec_J1J3;
    BCAO_6_3T = BCAO_meas_spec_J1J3;
    BCAO_6_4T = BCAO_meas_spec_J1J3;
    BCAO_6_5T = BCAO_meas_spec_J1J3;
    BCAO_27_3T = BCAO_meas_spec_J1J3;

    BCAO_6_2T.fname=f_ei6_2T;
    BCAO_6_3T.fname= f_ei6_3T;
    BCAO_6_4T.fname= f_ei6_4T;
    BCAO_6_5T.fname= f_ei6_5T;
    BCAO_27_3T.fname=f_ei27_3T;

    BCAO_6_2T.load_data();
    BCAO_6_3T.load_data();
    BCAO_6_4T.load_data();
    BCAO_6_5T.load_data();
    BCAO_27_3T.load_data();


    %Define masking for each Ei
    %No futher masking needed for Ei=6

    %Ei=27 needs phonon masking. 
    bad_i_27 = (BCAO_27_3T.E>15 & BCAO_27_3T.Q>-0.8 & BCAO_27_3T.Q<0.8 | BCAO_27_3T.E<1.5);
    bad_i_27_2 = (BCAO_27_3T.E>16 & abs(BCAO_27_3T.Q)>1.2);
    bad_i_27_3 = (abs(BCAO_27_3T.Q)>1.7);
    bad_i_27 = (bad_i_27 | bad_i_27_2 | bad_i_27_3);
    BCAO_27_3T.mask_meas_ind(bad_i_27);

    %Rather than fitting for a linear background later, it is easier to remove
    %a linear background using a simple average. 
    measurements = [BCAO_6_2T BCAO_6_3T BCAO_6_4T BCAO_6_5T BCAO_27_3T];
    for i=1:length(measurements)
        dataset = measurements(i);
        if dataset.Ei==6
            elim=[0.6,1.5];
            qlim=[-0.9,0.4];
        elseif dataset.Ei==27
            elim=[4,10];
            qlim=[-1.8,-1];
        end
        %Comment this line to skip this process.
        dataset.remove_simplebkg(qlim,elim);
    end

    %Put exchange params here, the initialize the structure
    %Or input params manually, below should be from a previous fit.
    measurements = [BCAO_6_3T BCAO_27_3T];
    measurements_field=measurements;
    J1xfit = J1x;
    J1zfit = J1z;
    J3xfit = J3x;
    J3zfit = J3z;
    gab_fit=gab;
    

    for i=1:length(measurements_field)
        dataset=measurements_field(i);
        dataset.initialize_exchange(J1xfit,J1zfit,0,0,J3xfit,J3zfit,0,0,0,0,0,0,0,0,3,gab_fit,gab_fit,2.0);
        dataset.generate_struct();
        dataset.calc_spec();
        dataset.calc_chisqr();
        fprintf('Chisqr = %d',dataset.chisqr)
    end

    % Plot for both Ei=6 and Ei=27 configs

    figure;
    spec6 = calc_BCAO_spec(BCAO_6_3T.structure,{[-1 0 3] [1 0 3] 200},linspace(0,6,300),BCAO_6_3T.res_mat);
    sw_plotspec(spec6,'legend',false);
    cbar = colorbar();
    ylabel(cbar,'I (a.u.)','FontSize',16);
    ylim([0 6]);
    xticklabels({'-1','-3/4','-1/2','-1/4','0','1/4','1/2','3/4','1'})
    colormap(viridis)
    ylabel('$\hbar\omega$ (meV)','Interpreter','latex','FontSize',16);
    xlabel('(H00) (r.l.u.)','FontSize',16);
    titlestr = sprintf('J=%.2f K=%.2f G=%.2f Gp=%.2f',J1xfit,J1zfit,J3x,J3z);
    title(titlestr);
    caxis([0 10]);

    figure;
    spec27 = calc_BCAO_spec(BCAO_27_3T.structure,{[-2 0 3] [2 0 3] 200},linspace(0,25,300),BCAO_27_3T.res_mat);
    sw_plotspec(spec27,'legend',false);
    cbar = colorbar();
    ylabel(cbar,'I (a.u.)','FontSize',16);
    ylim([0 20]);
    xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'})
    colormap(viridis)
    ylabel('$\hbar\omega$ (meV)','Interpreter','latex','FontSize',16);
    xlabel('(H00) (r.l.u.)','FontSize',16);
    titlestr = sprintf('J=%.2f K=%.2f G=%.2f Gp=%.2f',J1xfit,J1zfit,J3x,J3z);
    title(titlestr);
    caxis([0 5]);


    % Overplot the line dispersion onto the measurement
    figure;
    scatter(BCAO_6_3T.Q,BCAO_6_3T.E,55,BCAO_6_3T.I,'filled','s');
    sw_plotspec(spec27,'legend',false,'mode','disp','x0',[-2,-1],'color','w');
    ylim([0 6])
    xlim([-1 1])
    colormap(viridis)
    caxis([0 1.5e-3])
    titlestr = sprintf('J=%.2f K=%.2f G=%.2f Gp=%.2f',J1xfit,J1zfit,J3x,J3z);
    title(titlestr);%BCAO_27_3T.plot_measurement()


    figure;
    scatter(BCAO_27_3T.Q,BCAO_27_3T.E,55,BCAO_27_3T.I,'filled','s');
    sw_plotspec(spec27,'legend',false,'mode','disp','x0',[-2,-1],'color','w');
    ylim([0 20])
    xlim([-2 2])
    colormap(viridis)
    caxis([0 1.5e-4])
    titlestr = sprintf('J^1_{xy}=%.2f J^1_z=%.2f J^3_{xy}=%.2f J^3_z=%.2f',J1xfit,J1zfit,J3x,J3z);
    title(titlestr);%BCAO_27_3T.plot_measurement()
end
