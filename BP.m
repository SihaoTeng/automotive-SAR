% Back projection function

function data = BP(data)
   data.image = zeros(size(data.x_mat));
   
   for N = 1:data.Nch
       sRC = data.sRC(:,:,N);
        image = zeros(size(data.x_mat));
        
        % phase compensation matrix
        f = zeros(1,size(sRC,2));
        for fastt = 1:size(sRC,2)
            f(fastt) = data.fc + (fastt-1)*data.B/(size(sRC,2)-1);
        end
        phcomp = zeros(size(sRC,1),size(sRC,2));
        for row = 1:size(sRC,1)
            for col = 1:size(sRC,2)
                phcomp(row,col)=exp(1j*2*pi*f(col)*2*data.tardis(row)*cos(data.azi)/data.c);
            end
        end
        % Operating phase compensation 
%         sRC = sRC.*phcomp;

        Nfft = size(sRC,2);
%         r_vec = linspace(data.Rmin(N),data.Rmax(N),Nfft);
    
        for tau = 1:data.N_pulse    % slow time index
            r_vec = linspace(data.Rmin(tau,N),data.Rmax(tau,N),Nfft);
            rc = sRC(tau,:);
            dR = sqrt((data.Antx(N,tau)-data.x_mat).^2+(data.Anty(N,tau)-data.y_mat).^2);
            I = find(and(dR > min(r_vec),dR < max(r_vec)));
    
            image(I) = image(I) + interp1(r_vec,rc,dR(I),'linear').*exp(1j*4*pi*data.fc*dR(I)/data.c);
        end
     data.image =  data.image + image;
   end
return

