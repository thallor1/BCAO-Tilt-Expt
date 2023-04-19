%First import data
close('all')
pool_count=8; %number of threads for parallelization

%Best fit params from PNAS paper
J1x=-7.6;J1z=-1.2;J1px=0;J1pz=0;J3x=2.5;
J3z=-0.85;D=0.1;Ep=-0.1;F=0;G=0;D3=0;E3=0;F3=0;
G3=0;gx=5.0;gy=5.0;gz=2.0;

H = 7.0;


BCAO = gen_struct_BCAO_xxz(J1x,J1z,J1px,J1pz,J3x,J3z,D,Ep,F,G,D3,E3,F3,G3,H,gx,gy,gz);
%{
nQ = 201;
nE = 501;
Qhv = linspace(-1.5,1.5,nQ);
Qkv =  0;
Qlv =  linspace(-6.5,6.5,nQ);


[Qh,Qk,Ql]=ndgrid(Qhv,Qkv,Qlv);
Q = [Qh(:) Qk(:) Ql(:)]';

spec = BCAO.spinwave(Q);

Ecut = [3.25,3.75];

Ev = linspace(0,7,nE);

spec = sw_egrid(spec,'component','Sxx+Syy+Szz','Evect',Ev);
spec = sw_instrument(spec,'dE',0.07);
spec3D = reshape(spec.swConv,nE-1,nQ,nQ);



Eidx = find(Ev>Ecut(1) & Ev<Ecut(2));
figure;
cut1 = squeeze(sum(spec3D(Eidx,:,:),1))/numel(Eidx)/(Ev(2)-Ev(1));
imagesc(Qhv,Qlv,cut1);
set(gca,'Ydir','normal')
xlabel('(00L) (r.l.u.)')
ylabel('(H00) (r.l.u.)')
%caxis([0 3])
colorbar
%}
%%
%Replicate the spaghetti plot from expt, replacing the -HH0 direction with
%00L
res_file = 'Data Files/hys_res_ei6_fc300.txt';
res_table = readtable(res_file);
res_matrix=res_table{:,:};
res_matrix(:,2)=res_matrix(:,2).*3;

ratio = 5.417;
nQ =300;
%Original orientation below
spec = BCAO.spinwave({[1/2 0 0] [1/3 1/3 0] [0 1 0] [0 1/2 0] [1/3 1/3 0] [1/2 1/2 0] [1 0 0] [1/2 0 0] nQ});
spec = sw_neutron(spec);
spec = sw_egrid(spec,'component','Sperp','Evect',linspace(0,6,201));
spec = sw_instrument(spec,'dE',res_matrix,'norm',false,'fid',0,'dQ',0.02);

avgSpec = zeros(size(spec.swConv));
%Take a random sample of sample orientations such that 00L is always in the
%scattering plane. 
nRot = 101;
delTheta = rand(nRot,1)*2.0*3.1415;

%%
%Rotated frame read in from files

flist = dir('Path_indices_hkl_rotated/*.csv');
%Example file 
if isempty(gcp('nocreate'))
    parpool(pool_count);
end
parfor i=1:length(flist)
    f = flist(i);
    fname = f.name;
    fname_split = split(fname,'_');
    phi = str2double(fname_split(3));
    phi = pi*phi/180.0;
    %Phi = 0 corresponds to zero in-plane rotation, field is along the 010
    %direction
    Hinit = [0;H*sin(85.0*3.1415/180.0);H*cos(85.0*3.1415/180.0)];

    %Hinit = [0;H*sin(78.0*3.1415/180.0);H*cos(78.0*3.1415/180.0)];
    Rz = [cos(-phi) -sin(-phi) 0; sin(-phi) cos(-phi) 0; 0 0 1];
    Hrot = Rz*Hinit;
    Hx = Hrot(1);
    Hy = Hrot(2);
    Hz = Hrot(3);

    qpath_f = readtable(strcat('Path_indices_hkl_rotated/',fname));


    BCAO_rot=gen_struct_BCAO_xxz_vecH(J1x,J1z,J1px,J1pz,J3x,J3z,D,Ep,F,G,D3,E3,F3,G3,Hx,Hy,Hz,gx,gy,gz);

    Qpathrot = {table2array(qpath_f(1,:)) table2array(qpath_f(2,:)) ...
        table2array(qpath_f(3,:)) table2array(qpath_f(4,:)) table2array(qpath_f(5,:)) ...
        table2array(qpath_f(6,:)) table2array(qpath_f(7,:)) table2array(qpath_f(8,:)) nQ};
    specrot = BCAO_rot.spinwave(Qpathrot);

    specrot = sw_neutron(specrot);
    specrot = sw_egrid(specrot,'component','Sperp','Evect',linspace(0,6,201));
    specrot = sw_instrument(specrot,'dE',res_matrix,'norm',false,'fid',0,'dQ',0.02);

    avgSpec = avgSpec + specrot.swConv/length(flist);
    disp(i)
end
%%
spec.swConv=avgSpec;
figure;
sw_plotspec(spec,'mode','color','dE',0.1);
ylim([0 7])
caxis([0.0 1])

%Save to a csv file. 
writematrix(spec.swConv,'Tilted_spectra_files/BCAO_j1j3_swConv_tilt.csv')
writematrix(spec.hkl,'Tilted_spectra_files/BCAO_j1j3_hkl_tilt.csv')
writematrix(spec.Evect,'Tilted_spectra_files/BCAO_j1j3_omega_tilt.csv')