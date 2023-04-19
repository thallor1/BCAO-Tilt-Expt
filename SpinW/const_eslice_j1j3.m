%First import data
close('all')
pool_count=14; %number of threads for parallelization

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
measurements = [BCAO_6_3T BCAO_27_3T];
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
%After removing background, remove lower band for the Ei27 measurement
bad_i_27 = (BCAO_27_3T.E<8);
BCAO_27_3T.mask_meas_ind(bad_i_27);

%% Generate exchange 

%Below are testing params
J1x=-7.0;
J1z=-2.0;
J3x=-5-J1x;
J3z=-2-J1z;

%Below are the parameters in the draft
J1x=-7.7;
J1z=-1.20;
J3x=2.6;
J3z=-0.85;

D=0.1;
Ep=-0.1;
gx=5.0;
gy=5.0;
gz=2.0;

BCAO_6_3T.initialize_exchange(J1x,J1z,0,0,J3x,J3z,D,Ep,0,0,0,0,0,0,3,gx,gy,gz);
BCAO_6_3T.generate_struct();

bcao_struct = BCAO_6_3T.structure;

%%

nQ = 201;
nE = 501;
Qhv = linspace(-1,1,nQ);
Qkv = linspace(0,2,nQ);Qlv = 0;
[Qh, Qk, Ql] = ndgrid(Qhv,Qkv,Qlv);

% Create a list of Q point, with dimensions of [3 nQ^2].
Q = [Qh(:) Qk(:) Ql(:)];

spec = sq.spinwave(Q);

Ev = linspace(0,5,nE);
spec = sw_egrid(spec,'component','Sxx+Syy+Szz','Evect',Ev);
spec = sw_instrument(spec,'dE',0.1);