%This function reads a gmsh file and return three matrices.
%         nodes:   a matrix of nnodes x 3
%                  nodes(i,j) : coordenate j (either x,y or z) of node i
%         lines:   a matrix of nlines x 3  of boundary segments 
%                  lines(i,1)  index for boundary conditions
%                  lines(i,2)  initial segment node index
%                  lines(1,3)  final segment node index
%         elems:   a matrix of triangle conectivities
%                  elems(i,1)  node one of triangle i
%                  elems(i,2)  node two of triangle i
%                  elems(i,3)  node three of triangle i
%
% You have to create a triangulation and save the mesh using gmesh. Gmesh is a mesh generator program 
% that produces one, two and three dimensional unstructures meshes. FOr the boundary conditions create 
% Physical lines  and set a label. Do not set Physical points
% The official gmesh page where you can download the program is
%   http://www.geuz.org/gmsh/
% Tutorial for 2D meshes:
%   http://ffep.sourceforge.net/Download/gui_tutorial.pdf
%
% Manuel Garcia
% Eafit University. 2010

function [nodes, lines, elems]=readGmsh(archivo)


fid=fopen(archivo);

[fid,message]	= fopen(archivo,'r');
if(fid < 0)
    fprintf('Filename: %s. \n %s\n',archivo,message)
    
    error('error opening file:  ')
    exit (0)
end


%
%   Leer encabezado 
%
fscanf(fid,'$MeshFormat\n',1);
mline = fgetl(fid);
mline = fgetl(fid);
mline = fgetl(fid);

if strcmp(mline,'$PhysicalNames') 
  found =0 ;
  while  ~ found 
    mline = fgetl(fid);
    found = strcmp(mline,'$EndPhysicalNames');
  end
  mline = fgetl(fid);
end


if strcmp(mline,'$Nodes')   

   nnodes = fscanf(fid,'%d\n',1);
   nodes=zeros(nnodes,3);

   for i = 1:nnodes
       aux = fscanf(fid,'%d',1);
       nodes(i,:) = fscanf(fid,'%f %f %f\n',3 );
   end
   mline = fscanf(fid,'%s\n',1); %$EndNodes
  
end

mline = fscanf(fid,'%s\n',1); %$Elements

if strcmp(mline,'$Elements')   

   nelems = fscanf(fid,'%d\n',1);
   FullElems=zeros(nelems,5);
   nlines = 0;
   ntri = 0;
   
   for i = 1:nelems
       A = fscanf(fid,'%d %d %d',3 ); %  elemnumber elemtype nproperties
       
       props =fscanf(fid,'%d',A(3) ); %reads nproperties
       
       if A(2) == 1
           nlines = nlines +1 ;
           B= fscanf(fid,'%d %d\n',2 );
           FullElems(i,1) =A(2); % line
           FullElems(i,2) =props(1); %PhysicalProp
           FullElems(i,3) =B(1);
           FullElems(i,4) =B(2);
       else
           ntri = ntri+1;
           B= fscanf(fid,'%d %d %d\n',3 );
           FullElems(i,1) =A(2); %   triangle
           FullElems(i,2) =props(1); %PhysicalProp
           FullElems(i,3) =B(1);
           FullElems(i,4) =B(2);
           FullElems(i,5) =B(3);
       end
       FullElems(i,:) ;   
   end
   lines = zeros(nlines,3);
   elems = zeros(ntri,3);
   for i=1:nlines
       lines(i,:) = FullElems(i,2:4);
   end
   for i=1:ntri
       elems(i,:) = FullElems(i+nlines,3:5);
   end
  
end