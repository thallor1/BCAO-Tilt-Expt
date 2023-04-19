function chisqr_tot = total_chisqr(bcao_measurements,J1x,J1z,J3x,J3z,D,E,H,gx,gy,gz)
    % Provided intialized BCAO measurements (lidx and data already
    % generated) returns a total chisqr for all of the measurements. Also
    % requires a set of exchange parameters.
    chisqrs = zeros(length(bcao_measurements),1);
    for i=1:length(bcao_measurements)
        data = bcao_measurements(i);
        data.initialize_exchange(J1x,J1z,0,0,J3x,J3z,D,E,0,0,0,0,0,0,H,gx,gy,gz);
        data.generate_struct();
        try
            data.calc_spec();
            c2 = data.chisqr;
            %Weigh the Ei=27 meV measurement to be the same as the sum of the
            %Ei=6 measurements.
            totcount = length(bcao_measurements);
            chisqrs(i)=c2;
        catch
            chisqrs(i)=1e4;
        end
    end
    chisqr_tot = sum(chisqrs,'all')/length(bcao_measurements);
end