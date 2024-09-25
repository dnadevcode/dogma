A = [ 1 400 1 1;...
    100 500 2 0.95];

B = [5000 5400 2 0.97;...
    5100 5500 1 0.94];


import Core.update_sf_barset;
[Aupd] = update_sf_barset(A, B(1,4)/A(1,4), 2);

selBarId = 1;
Aupd(:,1:2) = Aupd(:,1:2)-Aupd(selBarId,1)+B(selBarId,1); %
