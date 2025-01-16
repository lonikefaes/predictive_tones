function FV = load_mesh(file,info)

mat          =  info.mat;
mesh         =  read_vtkMesh(file);
points       = (inv(mat)*[mesh.points';ones(1,size(mesh.points,1))])';
FV.vertices  =  points(:,1:3);
FV.faces     =  mesh.triangles;