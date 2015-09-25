## this script fits a covariance curve to the data set

filename = "H:/Global Analyser/TP-G1-01 Kriging_R_gstats/Data/brady02.csv";
M = csvread(filename);

## subtrack first time 
M(:,1)=M(:,1)-M(2,1);

data = M(2:end,2:end);
## subtract first settlement values 
data = data - data(1,:);
first =  data(90,:);

## load location 
locfile = "C:/dev/pyCode/Passer/locs.csv";
N = csvread(locfile);
sensor_locs = N(2:end,2:3)';
n_sensor_locs = size(sensor_locs,2);
dists= zeros(n_sensor_locs);

## find distance between each pair of sensors
for i = 1:n_sensor_locs
  for j=1:i
    dists(i,j) = dist([sensor_locs(1,i),sensor_locs(2,i)],[sensor_locs(1,j),sensor_locs(2,j)]);
    dists(j,i) = dists(i,j);
   end
end

## places distances in matrix into incrementing intervals of size - box_size in new matrix - boxes

box_size = 10;
maxunit = ceil(max(max(dists))/box_size)*box_size;
boxes = ceil(dists/box_size);
table = 0:box_size:maxunit;
means = zeros(1,size(table,2)-1);
count = zeros(1,size(table,2)-1);
vari = zeros(1,size(table,2)-1);



for i = 1:(maxunit/box_size)
  ## number of distance pairs in each interval 
  count(i) = nnz(boxes==i);
  ## mean distance value for each interval 
  means(i) = mean(dists(boxes==i));
  ## row and column indexes for each distance pair in boxes 
  [row col] = find(boxes==i);
  ## works out variance estimate
  vari(i) = sum((first(row)-first(col)).^2)/size(row,1);
end

## Mean Squared error function to minimize 
MSE = (@(b) (vari - b(1)-b(1)*b(3)+gauCov2(b(1),b(2),b(3),means)).^2*count');
## searches for minimum
xa = fminsearch(MSE, [0.05  50 0.3]);


## plots the solution to check it looks okay. 
t = 1:1:100;
plot(t,xa(1)*(1+xa(3))-gauCov2(xa(1),xa(2),xa(3),t));
hold on;
plot(means, vari, 'o')

xa = fminsearch(MSE, xa);
hold off;


time = M(2:end,1);
n_times = length(time);
time_diff = zeros(n_times);
time_diff= abs(repmat(time,1,n_times)-repmat(time',n_times,1));
[x y] = find(time_diff>800);

max_diff = 10;
box_size = 7;
maxunit = ceil(max(max(time_diff))/box_size)*box_size;
boxes = ceil(time_diff/box_size);
table = 0:box_size:maxunit;


