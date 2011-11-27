% Nacist soubor s namerenymi hodnotami
load pos.mat

% Ulozit jednotlive namerene vektory ve formatu, ktery
% umime nacist do Maple prikazem ImportMatrix. Napriklad
% pro vektor t0:
%     t0 := ImportMatrix("t0.mat", source=Matlab):
save('P0.mat', 'P0', '-ascii');
save('P1.mat', 'P1', '-ascii');
save('P2.mat', 'P2', '-ascii');
save('P3.mat', 'P3', '-ascii');
save('t0.mat', 't0', '-ascii');
save('t1.mat', 't1', '-ascii');
save('t2.mat', 't2', '-ascii');
save('t3.mat', 't3', '-ascii');