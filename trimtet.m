% 2025: 10/6 -- TRIMS a tetrahedral partition
% Uses (xr,yr,zr) as in points to decide which tets to trim


function  [x,y,z,TET,tetout] = trimtet(x,y,z,TET,xr,yr,zr);

nt = length(TET);
[it,B] = tsearchn([x,y,z],TET,[xr,yr,zr]);
tetout = setdiff([1:nt]',it);
tetin = setdiff([1:nt]',tetout);
TET = TET(tetin,:);

np=length(x);
rpts=[];
for i=1:np
    if sum(sum(ismember(TET,i)))==0
        rpts=[rpts;i];
    end
end

pts=[x y z];
pts(rpts,:)=[];

for i=1:length(TET)
    for j=1:4
        TET(i,j)=find(ismember(pts, [x(TET(i,j)),y(TET(i,j)),z(TET(i,j))],'rows'));
    end
end

x=pts(:,1); y=pts(:,2); z=pts(:,3);


