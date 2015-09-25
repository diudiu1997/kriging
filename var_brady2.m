## run this after find_cov

minGx = 85120;
maxGx = 85190;

minGy = 36500;
maxGy = 36700;

sensor_data = zeros(1,length(sensor_locs));

## makes mesh grid x and y

[gridx gridy] = meshgrid(minGx:2:maxGx , minGy:2:maxGy);
kriging_locs = [gridx(:)'; gridy(:)'];

## creates gauCov as an function variable with parameters found in find cov
funCov = (@(b) gauCov2(xa(1),xa(2),xa(3),b))

## does the kriging
[est, var] = krig_fun3(first, sensor_locs,kriging_locs,funCov);

## plots it
imagesc(reshape(var, size(gridy)))

hold on;


plot((sensor_locs(1,:)/2-minGx/2)+1,(sensor_locs(2,:)/2-minGy/2),'ro' )
imagesc(reshape(var, size(gridy)))
