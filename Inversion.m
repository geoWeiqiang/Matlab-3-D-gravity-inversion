
%% Step 00 
% clear,clc
clear,clc
QH8B='D:\program\Matlab2021a\bin\my\PROCODE\';
addpath([QH8B,filesep,'Gravity_Inversion_v1',filesep,'code']);
addpath(genpath([QH8B,filesep,'Utility']));
% 
%% Step 01 

mX = 0:5:200;
mY = 0:5:200;
mZ = 0:5:100;
Den = meshgrid(mX,mY,mZ);

Den=Den*0;

Den(16:24,10:18,9:13)=2;
Den(16:24,24:32,9:13)=3;
[hf,ha]=M3DSliceView(Den,mX,mY,mZ,''); 
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');

[Mindex,Mcoor,delta,dz]=M3D_meshCoder(mX,mY,mZ,[]);
Model=Den(:);

Slice_x=Den(20,:,:);  
Slice_x=reshape(Slice_x(:),length(mX),length(mZ));
figure;imagesc(mX,mZ,Slice_x')
xlabel('X/m')
ylabel('Z/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice X-Z');
axis equal

Slice_y=Den(:,30,:);  
Slice_y=reshape(Slice_y(:),length(mY),length(mZ));
figure;imagesc(mY,mZ,Slice_y')
xlabel('Y/m')
ylabel('Z/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice Y-Z');
axis equal 

Slice_z=Den(:,:,10);  
Slice_z=reshape(Slice_z(:),length(mY),length(mX));
figure;imagesc(mX,mY,Slice_z)
xlabel('X/m')
ylabel('Y/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice X-Y');
axis equal

%% Step 02 
dX = 0:5:200;
dY = 0:5:200;
dZ = 0;
[Dindex,Dcoor,~,~]=M3D_meshCoder(dX,dY,dZ,[]);

%% Step 03 
% inversion_kernelFile=[outPath,filesep,'Kernel'];  
% 默认该文件名，软件内部保存到观测数据文件夹中
IsLoadKernel=1;
inversion_kernelFile='kernel.mat';
if(IsLoadKernel~=0)
    Gg=B01_ComputeKernel(Dcoor,Mcoor,delta);
    disp(['核矩阵正在写入文件。。。。。',inversion_kernelFile])
    save (inversion_kernelFile, 'Gg','Mcoor','Dcoor','mX','mY','mZ');
    disp(['核矩阵系数已经存入文件：',inversion_kernelFile])
    %     clear('Gg','Z');
else
    disp('正在加载核矩阵....')
    load(inversion_kernelFile);
    disp('核矩阵加载完成！')
end

%% Step 04 
Data=Gg*Model;
% Data=awgn(Data,20);  
Data=double(Data);

save('Forward.mat','Data','Dcoor','Mcoor')

%% Step 05 
figure;
plot(Data);
Data2D=reshape(Data,length(dY),length(dX));
figure;
contourf(dX,dY,Data2D)
xlabel('X/m')
ylabel('Y/m')
colorbar;
title('Result of Forward + 2%高斯噪声');
axis equal

%% Step 06 
Nm=length(Model);
Par.Ztype=1;
Par.ebeta=2.5;  
Z=C01_DepthWeighted_kernel(Par.ebeta,Gg);
Par.L=[1,1,1];
[Wm,Wd]=C00_Weights(Data,Mindex,Nm,delta,Z,Par.L);
Par.Beta=600;  
Par.Threshold=1.0E-5;
Par.Iterate=400;
Par_bk=Par;
ic=0;
Mod_ref=Model*0.1;  
Mod_ini=0*ones(Nm,1); 
[Mod,phid,phim]=D01_inversion_unstr(Data,Gg,Wm,Wd,Par,Mod_ini,Mod_ref);
Mod=double(Mod);

grid on;
grid minor;
pbaspect([2 1 1]);

M3D=M3D_meshDecode(Mindex,Mod+Mod_ref);   
[hf,ha]=M3DSliceView(M3D,mX,mY,mZ,num2str(ic,'C_%d'));

for ic=1:0
    %Mod1 =Mod;
    %Mod1((Mod1-2).^2<1)=
    Mod_ini=Mod*0.;
    Mod_ref=Mod;
    [Mod,phid,phim]=D01_inversion_unstr(Data,Gg,Wm,Wd,Par,Mod_ini,Mod_ref);
    Mod=double(Mod);
    M3D=M3D_meshDecode(Mindex,Mod+Mod_ref);  
    [hf,ha]=M3DSliceView(M3D,mX,mY,mZ,num2str(ic,'C_%d'));
end

Slice_x=M3D(20,:,:);  
Slice_x=reshape(Slice_x(:),length(mX),length(mZ));
figure;imagesc(mX,mZ,Slice_x')
xlabel('X/m')
ylabel('Z/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice X-Z');
axis equal

Slice_y=M3D(:,30,:);  
Slice_y=reshape(Slice_y(:),length(mY),length(mZ));
figure;imagesc(mY,mZ,Slice_y')
xlabel('Y/m')
ylabel('Z/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice Y-Z');
axis equal 

Slice_z=M3D(:,:,10);   
Slice_z=reshape(Slice_z(:),length(mY),length(mX));
figure;imagesc(mX,mY,Slice_z)
xlabel('X/m')
ylabel('Y/m')
c1=colorbar;
set(get(c1,'title'),'string','g/cm^3');
title('Slice X-Y');
axis equal


                       