function [BCAO]=gen_struct_BCAO(J,K,G,Gp,J2,H,gx,gy,gz)
    %Parametrized JKG model as theta phi with an additional Gp allowed
    BCAO = spinw;
    warning off;
    Lattice=[5.007 5.007 23.491];
    % Make the lattice 
    
    J1 = J.*eye(3);
    
    %Kitaev terms
    Kxx = (K.*diag([1 0 0]));
    Kyy = (K.*diag([0 1 0]));
    Kzz = (K.*diag([0 0 1]));
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
    if H>0.55
        %KQSL / Field polarized phse
        S=[b_dir(1) b_dir(1) b_dir(1) b_dir(1) b_dir(1) b_dir(1); ...
            b_dir(2) b_dir(2) b_dir(2) b_dir(2) b_dir(2) b_dir(2); ...
            b_dir(3) b_dir(3) b_dir(3) b_dir(3) b_dir(3) b_dir(3)];
        S=0.5.*S;
        BCAO.genmagstr('mode','direct','k',k,'n',n,'S',S,'nExt',[1 1 1]);
    else
        %Need to move into coordinate system used in heylion paper:
        % y-axis is b, x-axis is a* and c-axis is c

        k = [1/3 0 -4/3];
        %See Reagnault Heylion paper for this
        %Expressions for moment along crystallographic directions
        %6 magnetic atoms in unit cell
        num_extra_cells=3;

        S = zeros(3,6*num_extra_cells);
        S(1,:)=1/2;
        S(2,:)=sqrt(3)/2;
        %Define indices of sublattice 1 and sublattice 2
        sublat_1 = [1;3;4];
        sublat_2 = [5;6;2];
        pos_sub1 = BCAO.matom.r(:,sublat_1);
        pos_sub2 = BCAO.matom.r(:,sublat_2); 
        %Add more unit cells.
        R_mat_A = zeros(3,3*num_extra_cells);
        R_mat_B = zeros(3,3*num_extra_cells);

        for i=1:num_extra_cells
            cell_i = i*3 -3;
            %Already did first cell
            pos_sub1_new = pos_sub1 + [i-1;0;0];
            pos_sub2_new = pos_sub2 + [i-1;0;0];
            R_mat_A(:,cell_i+1) = pos_sub1_new(:,1);
            R_mat_A(:,cell_i+2) = pos_sub1_new(:,2);
            R_mat_A(:,cell_i+3) = pos_sub1_new(:,3);
            %s_ferri script expects all A, then all B
            R_mat_B(:,cell_i+1) = pos_sub2_new(:,1);
            R_mat_B(:,cell_i+2) = pos_sub2_new(:,2);
            R_mat_B(:,cell_i+3) = pos_sub2_new(:,3);

        end
        k = [1/3 0 -4/3];
        mb = 2.5; % See Heylion
        mc=0.2;
        x = {R_mat_A R_mat_B k mb sublat_1 sublat_2};
        %define the function
        [S,k,n] = s_ferri(S,x);
        k=[1/3 0 0];
        BCAO.genmagstr('mode','direct','S',S,'n',n,'nExt',[num_extra_cells 1 1],'k',k); 
        
    end
    % Magnetic hamiltonian - generate bonds and couplings
    BCAO.gencoupling('maxDistance',6,'fid',0); %NN Co couplings

    
    % Define transformation matrix R from Kitaev coordinate to lattice coordinate
    %Below is above matrix but rotated into the unit cell as defined in
    %spinw
    R = [sqrt(2/3) 0 -1/sqrt(3);
        -1/sqrt(6) sqrt(2)/2 -1/sqrt(3);
        -1/sqrt(6) -sqrt(2)/2 -1/sqrt(3)];
    R_inv = inv(R);
    
    BCAO.addmatrix('label','Kxx','value',R_inv*Kxx*R ,'color','r');
    BCAO.addmatrix('label','Kyy','value',R_inv*Kyy*R,'color','g');
    BCAO.addmatrix('label','Kzz','value',R_inv*Kzz*R,'color','b');
    %BCAO.addmatrix('label','K','value',K,'color','r');
    % Heisenberg terms
    BCAO.addmatrix('label','J1','value',J1,'color','orange');
    if abs(J2)>0
        BCAO.addmatrix('label','J2','value',J2,'color','y');
    end
    %Off diagonal Gamma, Gamma prime terms

    BCAO.addmatrix('label','Gx','value',R_inv*([0 Gp Gp;Gp 0 G; Gp G 0])*R,'color','r');
    BCAO.addmatrix('label','Gy','value',R_inv*([0 Gp G;Gp 0 Gp; G Gp 0])*R,'color','g');
    BCAO.addmatrix('label','Gz','value',R_inv*([0 G Gp;G 0 Gp; Gp Gp 0])*R,'color','b');
     
    % Additional off diagonal terms allowed by trigonal distortion G'
    % X bond: G:{Sy Sz, Sz Sy} G':{Sx Sy, Sx Sz, Sy Sx, Sz Sx}
    % Y bond: G:{Sz Sx, Sx Sz} G':{Sy Sz, Sy Sx, Sz Sy, Sx Sy}
    % Z bond: G:{Sx Sy Sy Sx} G':{Sz Sx, Sz Sy, Sx Sz, Sy Sz}
    % These are as defined in literature coordinates.
    % Transform from literature to ours : Z_lit = X, X_lit = Y, Y_lit = Z
    % Add field: 
    
    b_dir = b_dir/norm(b_dir);
    BCAO.addmatrix('label','g_1','value',diag([gx gy gz]));
    BCAO.addg('g_1');
    BCAO.field(H.*b_dir);
    %BCAO.single_ion.g=[gx gy gz];
    
    % Add couplings - Kitaev 
    BCAO.addcoupling('mat','Kxx','bond',1,'subidx',[4 6 5]);
    BCAO.addcoupling('mat','Kyy','bond',1,'subidx',[7 9 8]);
    BCAO.addcoupling('mat','Kzz','bond',1,'subidx',[1 3 2]);
    % Heisenberg
    BCAO.addcoupling('mat','J1','bond',1);
    if abs(J2)>0
        BCAO.addcoupling('mat','J2','bond',2);
    end
    % Gamma + G'
    BCAO.addcoupling('mat','Gx','bond',1,'subidx',[4 6 5]); % Off-diagonal Gyz for N-neighbors
    BCAO.addcoupling('mat','Gy','bond',1,'subidx',[7 9 8]); % Off-diagonal Gyz for N-neighbors
    BCAO.addcoupling('mat','Gz','bond',1,'subidx',[1 3 2]); % Off-diagonal Gyz for N-neighbors
    

    %plot(BCAO,'range',[2 0.5 1],'bondRadius1',0.2,'bondMode','line',...
    %'bondLineWidth','lin','bondLinewidth0',0.2,'atomLegend',false);

    %Make ground state and hamiltonian compatible
    if H>=0.55 
        BCAO.optmagsteep('nRun',300,'fid',0);
    elseif H==0.4
        BCAO.anneal('initT',10,'random',true,'verbosity',0)
    end
end