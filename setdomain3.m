% setdomain 2025: 7/20
% (1) Read STL or point cloud [xp,yp,zp] -- stlread, scalecoud, mesh2pc, pcread

function [xp,yp,zp] = setdomain3;
% Read domain
fname = input('input filename (.stl, .dat, .ply) ');

% Read STL
if fname(end-2:end) == 'stl'
    stlData = stlread(fname);
    points = stlData.Points; np = length(points);
    xo = points(:,1); yo = points(:,2); zo = points(:,3);
    [x y z] = scalecloud(xo,yo,zo);
    TRI = stlData.ConnectivityList; nt = length(TRI);
    figure; trisurf(TRI,xo,yo,zo); colormap(winter);
    figure; trisurf(TRI,x,y,z); colormap(winter);
    shading interp; axis off;  lightangle(45,45); axis equal;

    
    % create point cloud
    np = 50000; %input('input np ');
    tic
    pc=mesh2pc([x,y,z],TRI,NumPoints=np);
    pts=pc.Location; np = length(pts);
    xp = pts(:,1); yp = pts(:,2); zp = pts(:,3);
    toc
    figure; plot3(xp,yp,zp,'.r'); axis equal;

else

    % Read point cloud (dat)
    if fname(end-2:end) == 'dat'
        data = dlmread(fname);
        np = data(1,1);
        xo = data(2:np+1,1); yo = data(2:np+1,2); zo = data(2:np+1,3);

    % Read point cloud (ply)
    elseif fname(end-2:end) == 'ply'
        pc = pcread(fname);
        pts = double(pc.Location);
        xo = pts(:,1); yo = pts(:,2); zo = pts(:,3);
    end

    A = unique([xo,yo,zo],'rows'); 
    [xp yp zp] = scalecloud(A(:,1),A(:,2),A(:,3));
    np = length(xp);
end

%[k,dis] = knnsearch([xp,yp,zp],[xp,yp,zp],'K',2);
%dis1 = dis(:,2);
%[min(dis1),sum(dis1)/np,max(dis1)]
