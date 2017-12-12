
function Error_du = errorH1P1P1(t_pos,du,duExacta_x,duExacta_y,geo)

global Time
tt = Time.time_vec(t_pos);

Error_du = 0;

for j=1:geo.nElement
    posdux = 6*j-5:6*j-3;
    posduy = 6*j-2:6*j;
    
    % Coord (x;y) de cada uno de los vertices del tríangulo
    coord = geo.coordinate(geo.element(j,:),:)';
      
    nCuad = 10;
    [X,Y,Wx,Wy] = triquad(nCuad,coord');
    test= ([coord;ones(1,3)])\[X(:),Y(:),ones(nCuad^2,1)]';

    duxAprox = test'*du(posdux);
    duyAprox = test'*du(posduy);
    
    diff_x =(duxAprox-duExacta_x(X(:),Y(:),repmat(tt,nCuad^2,1))).^2;
    diff_y= (duyAprox-duExacta_y(X(:),Y(:),repmat(tt,nCuad^2,1))).^2;
    
    Error_du= Error_du+ Wx'*reshape(diff_x+diff_y,size(Wx,1),size(Wy,1))*Wy;
end
