clc;
t=ncread('air.mon.mean.nc','air');
lon=ncread("air.mon.mean.nc",'lon');
lat=ncread("air.mon.mean.nc",'lat');
time=ncread("air.mon.mean.nc",'time');
LonSize=size(lon);LatSize=size(lat);TimeSize=size(time);

%----将温度t转化为时间-空间二维数组-----
T=zeros(LonSize(1)*LatSize(1),TimeSize(1));
for i=1:TimeSize(1)
    T(:,i)=reshape(t(:,:,i),[LonSize(1)*LatSize(1),1]);
end

%计算每年的平均,共73年
T_year=zeros(LonSize(1)*LatSize(1),73);
for i=1:73
    T_year(:,i)=mean(T(:,12*i-11:12*i),2);
end
%----年平均的平均----
T_mean=mean(T_year,2);
T_mean_1948_1990=mean(T_year(:,1:43),2);
T_mean_1990_2021=mean(T_year(:,43:73),2);
%距平化每年的数据
T_1948_1990=T_year(:,1:43)-T_mean_1948_1990;
T_1990_2021=T_year(:,43:73)-T_mean_1990_2021;
T_year=T_year-T_mean;


%----协方差矩阵----
X=T_year*T_year';
X2=T_1948_1990*T_1948_1990';
X3=T_1990_2021*T_1990_2021';
%求方差和特征向量
[V,D]=eig(X3);
lambda=wrev(diag(D));
V=fliplr(V);
V3=V(:,1:3);
PC=V3'*T_year;

lambda=lambda/sum(lambda);
rsum=zeros(LonSize(1)*LatSize(1),1);
rsum(1)=lambda(1);
for i=2:LonSize(1)*LatSize(1)
    rsum(i)=rsum(i-1)+lambda(i);
end