function spec=calc_BCAO_spec(struct,Q,E,res_mat)

    %Calculates a spin-wave spectra provided a spinw structure, Q, E, and
    %incident energy values.
    spec = struct.spinwave(Q,'hermit',false,'gtensor',true);
    spec = sw_neutron(spec);
    spec = sw_egrid(spec,'component','Sperp','Evect',E,'imagChk',false,'binType','ebin');
    spec = sw_instrument(spec,'dE',res_mat,'norm',false,'fid',0,'dQ',0.02);
    return
end