% Signal generation function

function data = signalGen(data)
    c = data.c;
    fc = data.fc;
    B = data.B;

    PRI = data.PRI;
    VPC_pos0 = data.VPC_pos0;    % original VPC position x coordinate and y coordinate
        
    FoV_min = data.FoV_min;   % FoV closest point
    FoV_max = data.FoV_max;  % FoV longest point
    vego = data.vego; % only move in x direction
    N_pulse = data.N_pulse;
    Nch = data.Nch;

    Rmin1 = sqrt((FoV_min(1)-VPC_pos0(1))^2+(FoV_min(2)-VPC_pos0(2))^2);   %produce artifact
    Rmax1 = sqrt((FoV_max(1)-VPC_pos0(1))^2+(FoV_max(2)-VPC_pos0(2))^2);   %produce artifact
    dt = 1/(2*B);   % step size smaller than 1/Nyquist (fast time)
    Nfast = length(2*Rmin1/c:dt:2*Rmax1/c-dt);    % fast time total of 184个

    data.sRC = zeros(N_pulse,Nfast,Nch);   % sampling points >= 184, chose 200
    data.Antx = zeros(Nch,N_pulse);
    data.Anty = zeros(Nch,N_pulse);
    data.Rmin = zeros(N_pulse,Nch);
    data.Rmax = zeros(N_pulse,Nch);
    data.tardis = zeros(1,N_pulse);

    for N_index = 1:Nch
        VPC_pos_ant = VPC_pos0 + [0,data.dy*(N_index-1)];
        sRC_single = zeros(N_pulse,Nfast);
        for tau = 1:N_pulse  % slow time index
            VPC_pos_new = VPC_pos_ant + (tau-1)*vego*PRI;
            data.Antx(N_index,tau) = VPC_pos_new(1);
            data.Anty(N_index,tau) = VPC_pos_new(2);
            Rmin = sqrt((FoV_min(1)-VPC_pos_new(1))^2+(FoV_min(2)-VPC_pos_new(2))^2);   %produce artifact
            Rmax = sqrt((FoV_max(1)-VPC_pos_new(1))^2+(FoV_max(2)-VPC_pos_new(2))^2);   %produce artifact
%             Rmin = sqrt((FoV_min(1)-data.Antx(N_index,1))^2+(FoV_min(2)-data.Anty(N_index,1))^2);
%             Rmax = sqrt((FoV_max(1)-data.Antx(N_index,1))^2+(FoV_max(2)-data.Anty(N_index,1))^2);
    
            Tp = 2*Rmax/(c*B);  % pulse duration
            dt = 1/(2*B);   % step size smaller than 1/Nyquist (fast time)
%             t = 2*Rmin/c:dt:2*Rmax/c-dt;    % fast time total of 184个
            t = linspace(2*Rmin/c, 2*Rmax/c-dt, Nfast); % sampling points >= 184 (fast time samples)
            data.Tp = Tp;
            data.t = t;

            % define moving target
            % target origin position = [x,y,A] x,y:coordinates, A: reflectivity
            target = data.target;   % 3 targets position and reflectivity:0-1

            % target motion
            target(2,1:2) = target(2,1:2)+data.vtarget*(tau-1)*PRI;
            data.tardis(tau) = data.vtarget(1)*(tau-1)*PRI;
            dx = abs(data.Antx(1,1)-data.target(3,1));
            dy = abs(data.Anty(1,1)-data.target(3,2));
            data.azi = atan(dy/dx);

            for tar_index = 1:size(target,1) 
                Rtar = sqrt((target(tar_index,1)-VPC_pos_new(1))^2+(target(tar_index,2)-VPC_pos_new(2))^2);
                TD = 2*Rtar/c;
                sRC_single(tau,:) = sRC_single(tau,:)+target(tar_index,3)*Tp*sinc(B*(t-TD))*exp(-1j*2*pi*fc*TD);
            end
            data.Rmin(tau,N_index) = Rmin;
            data.Rmax(tau,N_index) = Rmax;
        end
        data.sRC(:,:,N_index) = sRC_single;
    end 
return