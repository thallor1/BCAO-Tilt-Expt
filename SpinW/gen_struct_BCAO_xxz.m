function [BCAO]=gen_struct_BCAO_xxz(J1x,J1z,J1px,J1pz,J3x,J3z,D,E,F,G,D3,E3,F3,G3,H,gx,gy,gz)
    %Parametrized J1J3 with DM model as theta phi with an additional Gp allowed
    BCAO = spinw;
    warning off;
    Lattice=[5.007 5.007 23.491];
    % Make the lattice 

    BCAO.genlattice('lat_const',Lattice,'angled',[90 90 120],'spgr',148,'fid',0);
    BCAO.addatom('label','Co','r',[0.0; 0.0;0.17014],'S',1/2,'color','Blue');
    BCAO.addatom('label','Ba','r',[0.0;0.0;0.0117],'S',0,'color','lightGray');
    BCAO.addatom('label','As','r',[0.3333;0.66667;0.02145],'S',0,'color','DarkCyan');
    BCAO.addatom('label',{'O1','O2'},'color',{'Red','Red'},'S',[0,0],'r',[
        0.333 0.0163;0.66667 0.3375; 0.02145 0.11476])


    %Magnetic structure
    BCAO.lattice;
    k=[0 0 0]; % for fully polarized state
    b_dir = [-1/2 sqrt(3)/2 0]; %Neutron direction
    %b_dir = [1 0 0]; %THz direction
    b_dir = b_dir/norm(b_dir);
    n=[0 0 1];
    %Now generate the magnetic structure based on the field.
    %AFM1->AFM2->Polarized
    %H<0.35, H<0.55, H>0.55
    S=[b_dir(1) b_dir(1) b_dir(1) b_dir(1) b_dir(1) b_dir(1); ...
        b_dir(2) b_dir(2) b_dir(2) b_dir(2) b_dir(2) b_dir(2); ...
        b_dir(3) b_dir(3) b_dir(3) b_dir(3) b_dir(3) b_dir(3)];
    S=0.5.*S;
    BCAO.genmagstr('mode','direct','k',k,'n',n,'S',S,'nExt',[1 1 1]);

    BCAO.gencoupling('maxDistance',12,'fid',0); %NN, NNN Co couplings

    %NN Heisenberg, in-plane
    BCAO.addmatrix('label','J1','value',([J1x 0 0;0 J1x 0;0 0 J1z]),'color','b');
    %%NN Heisenberg, out-plane
    BCAO.addmatrix('label','J1p','value',([J1px 0 0;0 J1px 0;0 0 J1pz]),'color','g');
    %NNN Heisenberg, in-plane
    BCAO.addmatrix('label','J3','value',([J3x 0 0;0 J3x 0; 0 0 J3z]),'color','r');
    %DM Interaction
    BCAO.addmatrix('label','DM1','value',([D E G;E -D F;G F 0]),'color','m');
    BCAO.addmatrix('label','DM3','value',([D3 E3 G3;E3 -D3 F3;G3 F3 0]),'color','m');


    BCAO.addcoupling('mat','J1','bond',1);
    BCAO.addcoupling('mat','J1p','bond',6);
    BCAO.addcoupling('mat','J3','bond',3);
    BCAO.addcoupling('mat','DM1','bond',1);

    b_dir = b_dir/norm(b_dir);
    BCAO.addmatrix('label','g_1','value',diag([gx gy gz]));
    BCAO.addg('g_1');
    BCAO.field(H.*b_dir);

    %plot(BCAO,'range',[2 2 1],'bondRadius1',0.2,'bondMode','line',...
    %    'bondLineWidth','lin','bondLinewidth0',0.2,'atomLegend',false,'atomMode','mag');
    %Make ground state and hamiltonian compatible
    BCAO.optmagsteep('fid',0);
end