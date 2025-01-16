function get_fcn2(range,fcn_files,dm,TR,cutoff)
%addpath(genpath('C:\Analysis\spm8_10 - Copy\'));
path = cd;
cd([path,'\FCN\']);


%load range;
T = 218;
R = 6;


mx0 = range.r0.mx; 
my0 = range.r0.my; 
mz0 = range.r0.mz; 
mx5 = range.r5.mx; 
my5 = range.r5.my; 
mz5 = range.r5.mz;

nx0 = length(mx0);
ny0 = length(my0);
nz0 = length(mz0);

nx5 = length(mx5);
ny5 = length(my5);
nz5 = length(mz5);

interp_factor     = 5;  %!!!!!
FD5d =  cell(1,6);
FDd =  cell(1,6);
[FD5d{:}] = deal(zeros(T,nx5*ny5*nz5));
[FDd{:}]  = deal(zeros(T,nx0*ny0*nz0));
for k =1:length(fcn_files),
    P = [];
    FD = [];
    % % % read fucntional data:
    P   = spm_vol(fcn_files{k});

    nx0 = P(1).dim(1); 
    ny0 = P(1).dim(2);
    nz0 = P(1).dim(3);
    
   [X0 Y0 Z0]    = meshgrid(0:length(my0)-1,0:length(mx0)-1,0:length(mz0)-1);

   [Xr4 Yr4 Zr4] = meshgrid(linspace(0,length(my0)-1,length(my5)),...
                      linspace(0,length(mx0)-1,length(mx5)),...
                      linspace(0,length(mz0)-1,length(mz5)));
 
    method = 'nearest'; 
    FD   = zeros(P(1).dim(1:3));

    for i=1:length(P),
        j = 0;
        for p=1:P(i).dim(3),
            j= j+1;
            FD(:,:,j) = spm_slice_vol(P(i),spm_matrix([0 0 p]),P(i).dim(1:2),0);
        end,
        FD5   = interp3(X0,Y0,Z0,FD(mx0,my0,mz0),Xr4,Yr4,Zr4,method);
        FD0   = FD(mx0,my0,mz0); % 
    %    save(['FD_run',num2str(k),'_',num2str(i)],'FD0','FD5');
        disp(i);
        FD5d{k}(i,:) = FD5(:)';
        FDd{k}(i,:)  = FD0(:)';
    end

end 



% high-pass filter data


design_matrix = load(dm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(design_matrix,1)
    epoch_start(i,1) = design_matrix(i,1) - 2;
    epoch_end(i,1) = design_matrix(i,1) + 7 + 14;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for runs = 1:R
    FDd{runs} = filt_tc(FDd{runs},TR,cutoff,0);
    FD5d{runs} = filt_tc(FD5d{runs},TR,cutoff,0);
    for i = 1:size(design_matrix,1)
        Data0(:,:,i,runs)     =  FDd{runs}(epoch_start(i,1):epoch_end(i,1),:);
        DataR(:,:,i,runs)     =  FD5d{runs}(epoch_start(i,1):epoch_end(i,1),:);
    end
end

Data0F = squeeze(mean(squeeze(mean(Data0(:,:,:,1:3),4)),3));
Data0S = squeeze(mean(squeeze(mean(Data0(:,:,:,4:6),4)),3));
DataRF = squeeze(mean(squeeze(mean(DataR(:,:,:,1:3),4)),3));
DataRS = squeeze(mean(squeeze(mean(DataR(:,:,:,4:6),4)),3));


nx0 = length(mx0);
ny0 = length(my0);
nz0 = length(mz0);

nx5 = length(mx5);
ny5 = length(my5);
nz5 = length(mz5);

Data0FT = zeros(nx0,ny0,nz0,size(Data0,1));
DataRFT = zeros(nx5,ny5,nz5,size(Data0,1));
Data0ST = zeros(nx0,ny0,nz0,size(Data0,1));
DataRST = zeros(nx5,ny5,nz5,size(Data0,1));

for i = 1:size(Data0,1)
    Data0FT(:,:,:,i) = reshape(Data0F(i,:),[nx0 ny0 nz0]);
    DataRFT(:,:,:,i) = reshape(DataRF(i,:),[nx5 ny5 nz5]);
    Data0ST(:,:,:,i) = reshape(Data0S(i,:),[nx0 ny0 nz0]);
    DataRST(:,:,:,i) = reshape(DataRS(i,:),[nx5 ny5 nz5]);
end;

%
save('FD_avg_raw','Data0FT','DataRFT','Data0ST','DataRST');

cd(path);