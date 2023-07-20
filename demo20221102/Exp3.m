%%
clear;
clc;
close all;

%%
%parameters setting;
fs = 2e9; %sampling frequency
fin = 1e8; % input frequncy
nCyl = 2;



ups = 100;
T = 2/fin;
x =0:1/fs:T-1/fs;
x = x';






n = (fs/fin)*nCyl;


n_bin = (n+1)*ups;
bits = 12;

%%
interval_t = 0.1;
num_t = ceil((nCyl*1/fin)/(interval_t/fs));
% t=0:0.01/fs:nCyl*1/fin; %time index
t = 0:interval_t/fs:num_t*interval_t/fs;

x = fin*t/nCyl;

%%
y = upsample(x,ups);
h = ones(ups,1);
z = filter(h,1,y);

in1 = min(z);
inp = max(z);
in = ((z-in1)/(inp-in1));
% plot(z)

%%
edges = 0:10:(2^bits);
edges = (edges)/(2^bits);

ideal = histogram(in,edges);
% figure('Name','ideal histogram');
% stem(ideal)
% title('ideal histogram');
