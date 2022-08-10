%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Post Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A lot still TBD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output & Plot designated output values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_n = dt*(0:nMax-1);
time_nph = time_n+0.5*dt;
%  Loop over the output quantities:
for iout = 1:numOutputQty
    
    if OutputType(iout) <= 8
        switch OutputType(iout)
            case 1
                titleStr = ['Ex Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ex (V/m)';
                filename = ['ExProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 2
                titleStr = ['Ey Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ey (V/m)';
                filename = ['EyProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 3
                titleStr = ['Ez Probe ' num2str(iout,'%d')];
                ylabelStr = 'Ez (V/m)';
                filename = ['EzProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 4
                titleStr = ['Hx Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hx (A/m)';
                filename = ['HxProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 5
                titleStr = ['Hy Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hy (A/m)';
                filename = ['HyProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 6
                titleStr = ['Hz Probe ' num2str(iout,'%d')];
                ylabelStr = 'Hz (A/m)';
                filename = ['HzProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
            case 7
                titleStr = ['Voltage Probe ' num2str(iout,'%d') ];
                ylabelStr = 'V (V)';
                filename = ['VProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_n';
            case 8
                titleStr = ['Current Probe ' num2str(iout,'%d') ];
                ylabelStr = 'I (A)';
                filename = ['IProbe_' num2str(iout,'%d') '.txt'];
                dataOut(:,1) = time_nph';
                
        end
        
        % plot if requested
        if outputPlot(iout) > 0
            figure;
            plot(dataOut(:,1),OutputValue(:,iout),'Linewidth',2);

            title(titleStr);
            xlabel('t (s)');
            ylabel(ylabelStr);
        end

        % write data to txt file:
        [fid, msg] = fopen(filename,'w');

        dataOut(:,2) = OutputValue(:,iout);
        if fid > 0
            for no = 1:nMax
                fprintf(fid,'%d %d \n',dataOut(no,1),dataOut(no,2));
            end
            fclose(fid);
        else
            msg
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Near to Far Field Post Processing and Graphing!!!~!!!!
%---Conversion from cartesian L,N compoents to Spherical components-----%%%
for f=1:length(frequencies) %%%CHECKED%%%
    for theta=1:length(theta_angles)
        for phi=1:length(phi_angles)
            T=theta_angles(theta);
            P=phi_angles(phi);
           
            
            L_theta(f,theta,phi)=(Lx(f,theta,phi)*cos(T)*cos(P))+(Ly(f,theta,phi)*cos(T)*sin(P))-(Lz(f,theta,phi)*sin(T));
            L_phi(f,theta,phi)=(-Lx(f,theta,phi)*sin(P))+(Ly(f,theta,phi)*cos(P));
            
            N_theta(f,theta,phi)=(Nx(f,theta,phi)*cos(T)*cos(P))+(Ny(f,theta,phi)*cos(T)*sin(P))-(Nz(f,theta,phi)*sin(T));
            N_phi(f,theta,phi)=(-Nx(f,theta,phi)*sin(P))+(Ny(f,theta,phi)*cos(P));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------PERFORM Eqs 8.51, 8.52,8.54,8.55----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r=1:length(radial_distances)%%%checked
    R=radial_distances(r);
    for f=1:length(frequencies)
        k=(2*pi*frequencies(f))*sqrt(eps0*mu0)
        for theta=1:length(theta_angles)
            for phi=1:length(phi_angles)
                

                Etheta(r,f,theta,phi)=((-1i*k)*(1/(4*pi*R)))*(exp(-1i*k*R))*((L_phi(f,theta,phi)) + eta0*(N_theta(f,theta,phi)))
                Ephi(r,f,theta,phi)=((1i*k)*(1/(4*pi*R)))*(exp(-1i*k*R))*((L_theta(f,theta,phi)) - eta0*(N_phi(f,theta,phi)));

                Htheta(r,f,theta,phi)=(-1/eta0) * Ephi(r,f,theta,phi);
                Hphi(r,f,theta,phi)  = (1/eta0) * Etheta(r,f,theta,phi);

            end
    
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------PERFORM ATTENA GAIN CALCULATIONS----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S11=0; %?
%V1=1;
%Z01=50; %?
%Pin=(1-(abs(S11))^2)*(((abs(V1))^2)/(2*Z01))
Pin=1;  %figure out what power is fed into attenna

%%Initialze Gain Arrays

G_theta=zeros(r,f,theta,phi);
G_phi=zeros(r,f,theta,phi);
G_total=zeros(r,f,theta,phi);

%%%%%%%%Calculate Gains
for r=1:length(radial_distances)
    R=radial_distances(r);
    for f=1:F_dim
        for theta=1:T_dim
            for phi=1:P_dim
                G_theta(r,f,theta,phi)=((4*pi*R^2)/(Pin))*((abs(Etheta(r,f,theta,phi)))^2)/(2*eta0);
                G_phi(r,f,theta,phi)=((4*pi*R^2)/(Pin))*((abs(Ephi(r,f,theta,phi)))^2)/(2*eta0);
                G_total(r,f,theta,phi)= G_theta(r,f,theta,phi) + G_phi(r,f,theta,phi);

            end
        end
    end
end

%ETT=squeeze(abs(Etheta(1,1,:,1)));


%polarplot(theta_angles,ETT);
%rmax=log(max(Gtheta));



Gtheta=squeeze(G_total(1,1,:,1));
Gphi=squeeze(G_total(1,1,1,:));


polarplot(theta_angles,log(Gtheta));
rmax=log(max(Gtheta));
text(0, rmax/2, 'dB', 'horiz', 'center', 'vert', 'top', 'rotation', 0);
%text(0, rmax*1.4, 'Theta', 'horiz', 'center', 'rotation', 0);
%title('Gain Vs Theta, R=10m f=300Mhz');





%polarplot(phi_angles,log(Gphi));
%rmax=log(max(Gphi));
%text(0, rmax/2, 'dB', 'horiz', 'center', 'vert', 'top', 'rotation', 0);
%text(0, rmax*1.4, 'phi', 'horiz', 'center', 'rotation', 0);
%title('Gain Vs Phi, R=10m f=300Mhz');

%subplot(2,1,2); 
%polarplot(phi_angles,-20*log(Gphi));



