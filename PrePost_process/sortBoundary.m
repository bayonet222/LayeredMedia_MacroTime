%*********************************************************************
%*                        Sort Boundary                              *
%*********************************************************************
%
%This function sort the edges in each boundary, and identificates de
%direction of each one
%
%***------------------------------------
%*** Inputs:  geometry, Gamma, indBoundary = [indDirichlet, indNewmann]
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [NewmanEdges,DirichletEdges] = sortBoundary(geo,Gamma,indBoundary)

edge2element = geo.edge2element;
indDirichlet = indBoundary(1);
indNewman    = indBoundary(2);
nodes2edge = geo.nodes2edge;

%% Inicialization

GammaN = Gamma(Gamma(:,1)==indNewman,:);
GammaD = Gamma(Gamma(:,1)==indDirichlet,:);

iniciales = edge2element(:,1);

%%
%*********************************************************************
%*                                                                   *
%*                         Newmann Boundary                          *
%*                                                                   *
%*********************************************************************

nEdgesNew    = length(GammaN(:,1));
NewmanEdges1 = zeros(nEdgesNew,4);
NewmanEdges2 = zeros(nEdgesNew,4);

for i=1:nEdgesNew
    
    %% Cheaking both directions of the boundary
    
    % Clockwise
    pos1      = find(iniciales==GammaN(i,2));
    finales1  = edge2element(pos1,2);
    elemento1 = edge2element(pos1(finales1==GammaN(i,3)),3);
    
    if ~isempty(elemento1)
        NewmanEdges1(i,:) = [nodes2edge(GammaN(i,2),GammaN(i,3)) GammaN(i,2:3) elemento1];
    end
    
    % Anticlockwise
    pos2      = find(iniciales==GammaN(i,3));
    finales2  = edge2element(pos2,2);
    elemento2 = edge2element(pos2(finales2==GammaN(i,2)),3);
    
    if ~isempty(elemento2)
        NewmanEdges2(i,:) = [nodes2edge(GammaN(i,2),GammaN(i,3)) GammaN(i,2:3) elemento2];
    end
end
if sum(sum(NewmanEdges1))==0
    NewmanEdges = sortrows(NewmanEdges2,1);
else
    NewmanEdges = sortrows(NewmanEdges1,1);
end


%%
%*********************************************************************
%*                                                                   *
%*                       Dirichlet Boundary                          *
%*                                                                   *
%*********************************************************************

nEdgesDir       = length(GammaD(:,1));
DirichletEdges1 = zeros(nEdgesDir,4);
DirichletEdges2 = zeros(nEdgesDir,4);

for i=1:nEdgesDir
%    for i=1:62
    %% Cheaking both directions of the boundary
    % Clockwise
    pos1      = find(iniciales==GammaD(i,2));
    finales1  = edge2element(pos1,2);
    elemento1 = edge2element(pos1(finales1==GammaD(i,3)),3);
    
    if ~isempty(elemento1)
        DirichletEdges1(i,:)= [nodes2edge(GammaD(i,2),GammaD(i,3)) GammaD(i,2:3) elemento1];
    end
    
    % Anti-clockwise
    pos2      = find(iniciales==GammaD(i,3));
    finales2  = edge2element(pos2,2);
    elemento2 = edge2element(pos2(finales2==GammaD(i,2)),3);
    
    if ~isempty(elemento2)
        DirichletEdges2(i,:)= [nodes2edge(GammaD(i,2),GammaD(i,3)) GammaD(i,2:3) elemento2];
    end
end
if sum(sum(DirichletEdges1))==0
    DirichletEdges = DirichletEdges2;
%     DirichletEdges = sortrows(DirichletEdges2,1);
else
    DirichletEdges = DirichletEdges1;
%     DirichletEdges = sortrows(DirichletEdges1,1);
end
