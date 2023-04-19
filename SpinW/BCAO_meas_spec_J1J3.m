classdef BCAO_meas_spec_J1J3 < handle
    %Object to store both the experimental and spinW information about a
    %particular configuration for BCAO in the HYS experiment 
    properties
        %Below values are from the experiment.
        fname %filename with measurement
        Q %1D list of measured Q (H00) in this experiment
        E %1D list of measured Energies
        I %1D list of measured intensities
        Err %1D list of measured errors
        res_mat=false %inelastic resolution
        Ei %incident energy of measurement 
        H % Field, also used in hamiltonian below
        label % label of this measurement. 
        
        %The below cannot be changed or the indexing breaks.
        Qcalc %Q List used to calculate the spin wave spectrum
        Ecalc %Energy list used to calculate the spin wave spectrum. 
        Qcalc_list %List of unique Q-values- different format than spinW notation
        lidx=false %Ultimately maps the flatted spinw spectra to the measurement. 
                %Should only be calculated once.
        
        %Below are the parameters relevant to the Hamiltonian and fitting.
        structure % SpinW structure object containing magnetic structure
        J1x %= -14.5;
        J1z %= -3.8;
        J1px%=0.0;
        J1pz%=0.0;
        J3x %= 7.2;
        J3z %= 2.3;
        D %= 1.5;
        Ep %= -1.7;
        F %= 0.0;
        G %= 0.0;
        D3% = 0.0;
        E3% = 0.0;
        F3% = 0.0;
        G3% = 0.0;
        gx % ~5
        gy % ~5
        gz % ~2 
        
        %Below are the calculated values
        spec % most recent spinw spectrum
        spec_scale % overall scaling factor for the calculated spectrum
        flatI % flattened spectra that matches the measured intensities. lidx is used for this
        flatQ % flattened Q
        flatE % flattened E
        chisqr % Chisqr between calculated spectra and measurement. 
    end
    methods
        function self = load_data(self)
            data = readtable(self.fname);
            if ~contains(self.fname,'Ei27') 
                self.Ei=6.0;
            else
                self.Ei=27.0;
            end
            self.I=table2array(data(:,1));
            self.Err=table2array(data(:,2));
            self.Q=table2array(data(:,3));
            self.E=table2array(data(:,4));
            %assign field based on file 
            if contains(self.fname,'2T')
                self.H=2.0;
            elseif contains(self.fname,'3T')
                self.H=3.0;
            elseif contains(self.fname,'=4T')
                self.H=4.0;
            elseif contains(self.fname,'5T')
                self.H=5.0;
            elseif contains(self.fname,'0p4T')
                self.H=0.4;
            elseif contains(self.fname,'=0T')
                self.H=0.0;
            end
            %Automatically remove bad values here. 
            good_i = (~isnan(self.I) & self.I>0 & self.Err>0 & ~isnan(self.Err));
            self.I = self.I(good_i);
            self.Err = self.Err(good_i);
            self.E = self.E(good_i);
            self.Q = self.Q(good_i);
            %Perform initial maksing of values based on Ei
            if self.Ei==6.0
                good_i = (self.E>0.35);
            elseif self.Ei==27.0
                good_i = (self.E>1.0);
            end
            self.I = self.I(good_i);
            self.Err = self.Err(good_i);
            self.E = self.E(good_i);
            self.Q = self.Q(good_i);
            %Initialize the calculated Q,E-values 
            self.Qcalc_list = unique(self.Q);
            self.Qcalc = {[min(self.Qcalc_list) 0 2] [0 0 2] [max(self.Qcalc_list) 0 2] length(self.Qcalc_list)};
            self.Ecalc = linspace(min(self.E),max(self.E),length(unique(self.E)));
        end
        
        function mask_meas_ind(self,bad_i)
            %Provided a set of indices to match, removes them from the
            %list. 
            good_i = (~bad_i);
            self.I = self.I(good_i);
            self.Err = self.Err(good_i);
            self.Q = self.Q(good_i);
            self.E = self.E(good_i);
        end
        
        function fig = plot_measurement(self)
            %Plots the stored measurement. 
            pointsize=180;
            scatter(self.Q,self.E,pointsize,self.I,'filled','s');
            titlestr = sprintf('BCAO Ei=%.1f H=%.1f T\n',self.Ei,self.H);
            title(titlestr)
            if self.Ei==6.0
                ylim([0,5.5])
            else
                ylim([0,20])
            end
            caxis([0 5.0*mean(self.I)])
            cbar = colorbar();
            ylabel(cbar,'I (a.u.)')
            xlabel('(H00)')
            ylabel('E (meV)')
        end
        
        function fig = plot_calc_spec(self)
            pointsize=180;
            scatter(self.flatQ,self.flatE,pointsize,self.flatI,'filled','s');
            disp(self.chisqr)
            titlestr=sprintf('BCAO Chi2=%.4f H=%.1f T\nJ1x=%.2f J1z=%.2f J3x=%.2f J3z=%.2f',self.chisqr,self.H,self.J1x,self.J1z,self.J3x, ... 
                self.J3z);
            title(titlestr)
            if self.Ei==6.0
                ylim([0,5.5])
            else
                ylim([0,20])
            end
            caxis([0 5.0*mean(self.I)])
            cbar = colorbar();
            ylabel(cbar,'I (a.u.)') 
            xlabel('(H00)')
            ylabel('E (meV)')
        end
        
        function remove_simplebkg(self,qlim_bkg,elim_bkg)
            box_i = (self.Q>qlim_bkg(1) & self.Q<qlim_bkg(2) & self.E>elim_bkg(1) & self.E<elim_bkg(2));
            I_bkg_all = self.I(box_i);
            Err_bkg_all = self.Err(box_i);
            weights = 1.0./Err_bkg_all;
            sum_weights = sum(weights,'all');
            mean_I_box = sum(I_bkg_all.*weights,'all')/sum_weights;
            self.I = self.I-mean_I_box;
        end
        function initialize_exchange(self,J1x,J1z,J1px,J1pz,J3x,J3z,D,Ep,F,G,D3,E3,F3,G3,H,gx,gy,gz)
            %Helper function to quickly assign relevant parameters.
            self.J1x=J1x; self.J1z=J1z; self.J1px=J1px; self.J1pz=J1pz; 
            self.J3x=J3x; self.J3z=J3z; self.D=D; self.Ep=Ep; self.F=F;
            self.G=G; self.D3=D3; self.E3=E3; self.F3=F3; self.G3=G3;
            self.gx=gx; self.gy=gy; self.gz=gz; self.H=H;
        end
        function generate_struct(self)
            %This function assumes that all exchange parameters have been
            %assigned. to the object. 
            BCAO = gen_struct_BCAO_xxz(self.J1x,self.J1z,self.J1px,self.J1pz,self.J3x,self.J3z,self.D,self.Ep, ... 
                self.F,self.G,self.D3,self.E3,self.F3,self.G3,self.H,self.gx,self.gy,self.gz);
            self.structure = BCAO;
        end
        
        function calc_spec(self)
            %Using the current exchange parameters, calculates LSWT
            %spectra. This is specifically for the first spectra. This
            %function will run the calculation, load the resolution matrix,
            % and then generate the lidx to map the flattened spectra onto
            % the 1D measured intensities. 
            if self.res_mat==false
                if self.Ei == 6.0
                    res_file = 'Data Files/hys_res_ei6_fc300.txt';
                elseif self.Ei==27
                    res_file = 'Data Files/hys_res_ei27_fc420.txt';
                else
                    disp('WARNING - only allowed Ei values are 6meV and 27 meV')
                end
                res_table = readtable(res_file);
                res_matrix=res_table{:,:};
                res_matrix(:,2)=res_matrix(:,2).*3;
                self.res_mat = res_matrix; %Will be used in future lswt
            end
            %Now that the resolution is handled, we can calculate a
            %spectrum
            lswt_spec = calc_BCAO_spec(self.structure,self.Qcalc,self.Ecalc,self.res_mat);
            self.spec = lswt_spec;
            %If the lidx indices have not been calculated, generate them
            %now. 
            if self.lidx==false
                h_calc = transpose(lswt_spec.hkl(1,:));
                e_calc = transpose(lswt_spec.Evect);
                e_calc = e_calc(1:length(e_calc)-1);
                e_calc = e_calc+abs(e_calc(2)-e_calc(1))/2.0;
                [Qc,Ec] = meshgrid(h_calc,e_calc);
                q_flat = reshape(Qc.',[],1);
                e_flat = reshape(Ec.',[],1);
                [~,lin_i]= min(abs(q_flat-self.Q.')+abs(e_flat-self.E.'));
                self.lidx=lin_i.'; %These will map the flattened spectra to the measurement.
                self.flatQ = q_flat(self.lidx);
                self.flatE = e_flat(self.lidx);
            end
            %Now that this is handled, store the flattened versions of the
            %calculated intensities and Q, E values into lists.
            specI = self.spec.swConv;
            %Scale the overall spectra to most closely match the measured
            %one 
            spec_flat = reshape(specI.',[],1);
            self.flatI = spec_flat(self.lidx);   
            scale_factor = get_theory_scaling(self.I,self.Err,self.flatI);
            self.flatI = self.flatI*scale_factor;
            self.spec_scale = scale_factor;
            self.spec.swConv = self.spec.swConv*self.spec_scale;
            %Now calculate a chisqr
            diff_sqr = abs(self.I - self.flatI).^2;
            errsqr = self.Err.^2;
            N = length(self.I);
            chisqr_spec = sum(diff_sqr./errsqr,'all')/N;
            if self.Ei==27
                chisqr_spec=chisqr_spec*3;
            end
            self.chisqr=chisqr_spec;
        end
        
        
        function calc_chisqr(self)
            % if a calculated spectra already exists, calculate the chisqr
            % between it and the measured spectra
            diff_sqr = abs(self.I - self.flatI).^2;
            errsqr = self.Err.^2;
            
            N = length(self.I);
            chisqr_spec = sum(diff_sqr./errsqr,'all')/N;
            if self.Ei==27
                chisqr_spec=chisqr_spec*3;
            end
            self.chisqr=chisqr_spec;

        end
    end
end
        