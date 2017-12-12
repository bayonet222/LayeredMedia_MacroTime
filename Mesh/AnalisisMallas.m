clc

%*********************************************************************
% Auxiliar code - MESH ANALYSIS
% Run after pre-process

figure('name', 'Micro - Scale')

%% Mesh
trimesh(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    zeros(size(Micro_geo.coordinate,1),1))
grid off
view(0,90)

hold on

[i,j,s] = find(Micro_geo.nodes2edge);
puestos = 1:Micro_geo.nnodes;
colors  = eye(3); %RGB

for e=1:length(s)%number of edges
    
    %% Number of vertex
    if puestos(i(e))~=0
        txt1 = sprintf('%i',i(e));
        text(Micro_geo.coordinate(i(e),1),Micro_geo.coordinate(i(e),2),txt1,'Color','b')
        puestos(i(e)) = 0;
    end
    if puestos(j(e))~=0
        txt2 = sprintf('%i',j(e));
        text(Micro_geo.coordinate(j(e),1),Micro_geo.coordinate(j(e),2),txt2,'Color','b')
        puestos(j(e)) = 0;
    end
    
    %% Number of edges
    p = [Micro_geo.coordinate(i(e),:)',Micro_geo.coordinate(j(e),:)'];
    pm  = (p(:,1)+p(:,2))./2;
    txt2 = sprintf('%i',s(e));
    text(pm(1),pm(2),txt2,'Color','b')
end
for e=1:Micro_geo.nElement %number of elements
    
    %% Number of element
    txt2 = sprintf('%i',e);
    pm = sum(Micro_geo.coordinate(Micro_geo.element(e,:),:))/3;
    text(pm(1),pm(2),txt2,'Color','k')
    
    coord = Micro_geo.coordinate(Micro_geo.element(e,:),:)';
    p  = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    normalvector = Micro_geo.normals{e};
    for n=1:3
        % EN CONTRA DEL RELOJ -- ASI DEBE SER SIEMPRE
        %normalvector(1:2,n) = (p(1:2,n)-p(3:4,n))'*[0,1;-1,0]/le(n);
        aaa = (p(1:2,n)+p(3:4,n))/2;
        q=quiver(aaa(1),aaa(2),0.05*normalvector(1,n),0.05*normalvector(2,n));
        
        q.Color = colors(n,:);
        q.LineWidth =0.75;
        q.MaxHeadSize = 1;
        
        %         % Numeration local vertex
        %         txt2 = sprintf(' %i-%i',e,n);
        %         xx= Micro_geo.coordinate(Micro_geo.element(e,n),1);
        %         yy= Micro_geo.coordinate(Micro_geo.element(e,n),2);
        %         text(xx+0.05*normalvector(1,n),yy+0.05*normalvector(2,n),txt2,'Color','k')
    end
    
end
set(gca,'xtick',[])
set(gca,'ytick',[])
