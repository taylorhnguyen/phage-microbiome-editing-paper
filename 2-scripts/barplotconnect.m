% things to pass this function:
% - data is the set of relative abundances, and should be ns x ng, where
%       ns is the number of samples and ng is the number of OTUs
% - time is a vector of the times of each sample
% - colors is a matrix of colors to use for the plot, will be generated
%       randomly if you input []
% - sc is the width of each bar separated by 1 day (default is 0.75 if you
%       input 0) 
% - trans is the transparency of the intermediate patches (default is 0.4
%       if you put -1)
% - last is the width of the last bar, which will be scaled by sc (default
%       is 1 if input is 0)

function [colors sc last] = barplotconnect(data,time,colors,sc,trans,last)
% table of data
%ns = 10; % number of samples
%ng = 20; % number of groups
%data = rand(ns,ng);
ns = size(data,1);
ng = size(data,2);
%time = [0 0.5 1:8];

if trans==-1
    trans = 0.4;
end

if sc==0;
    sc = 0.75;
end

if last==0
    last = 1;
end
   
% normalize the data
for i=1:ns
    data(i,:) = data(i,:)/sum(data(i,:));
end

% colors for each 
if numel(colors)==0
    colors = rand(ng,3);
end

% determine width of each
dt = diff(time);
t1 = time(1:end-1)+sc*dt;
%last = 1;
t1(end+1) = time(end)+sc*last;

%% draw patches
%figure;
hold on;
for i=1:ns
    tot = 0;
    for j=1:ng
        % make patch
        X = [time(i) t1(i) t1(i) time(i)];
        Y = [tot tot tot+data(i,j) tot+data(i,j)];
        patch(X,Y,colors(j,:));
        sum0(i,j) = tot;
        sum1(i,j) = tot+data(i,j);
        tot = tot+data(i,j);
    end
end

% now make patches
for i=1:ns-1
    for j=1:ng
        X = [t1(i) time(i+1) time(i+1) t1(i)];
        Y = [sum0(i,j) sum0(i+1,j) sum1(i+1,j) sum1(i,j)];
        p = patch(X,Y,colors(j,:));
        p.FaceAlpha = trans;
    end
end

xlim([time(1) t1(end)]);
ylim([0 1]);

    