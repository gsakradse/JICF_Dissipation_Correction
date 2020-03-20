function [mdx, mdv2] = PIV_Struc_FunY(Xmat, vel)

% Function computes the longitudinal structuere funciton in the x direction
% with the PIV data which is re-oriented.

% Making it loop through to accept columns of SCmat 2/28/20

% Inputs are Xvec and a single velocity snapshot.

% vel = Vr(:,:,50);
% x = Ymat(:,1)';



ninc = 200;
lag = 1:ninc;
tol = 0.15;


for kk = 1:size(Xmat,2)
    
    x = Xmat(:,kk);
    
    mdx = zeros(ninc,length(vel(1,:)));
    mdv2 = mdx;

    for jj = 1:length(vel(1,:))

        v = vel(:,jj)';

        for ii = lag

            % Find increments
           dx = x(1+ii:end) - x(1:end-ii);
           dv = v(1+ii:end) - v(1:end-ii);

           %mean increment
           meandx = mean(dx);

           ind = find(abs(dx-meandx)<tol);

           mdx(ii,jj) = mean(dx(ind));

           mdx(ii,jj) = nanmean(dx(ind))';
           mdv2(ii,jj) = nanmean(dv(ind).^2)';

        end   
    end
end