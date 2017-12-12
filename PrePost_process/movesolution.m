function Micro_Sol = movesolution(X,Micro_geo,Micro_solution1,Micro_solution2,Micro_Sol)

n = size(Micro_geo.coordinate,1);

solZ1 = Micro_solution1(repmat(X(1),n,1),repmat(X(2),n,1),...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));

solZ2 = Micro_solution2(repmat(X(1),n,1),repmat(X(2),n,1),...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));

diff1 = max(solZ1)- max(Micro_Sol.Pres1_cont);
Micro_Sol.Pres1_cont = Micro_Sol.Pres1_cont+diff1;

diff2 = max(solZ2)- max(Micro_Sol.Pres2_cont);
Micro_Sol.Pres2_cont = Micro_Sol.Pres2_cont+diff2;