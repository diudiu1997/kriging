## this is the main kriging function
## Takes a functional input covFun

function [est,var] = krig_fun3(sensor_data, sensor_locs,kriging_locs,covFun)

  n_sensor_locs = size(sensor_locs,2);
  n_kriging_locs = size(kriging_locs,2);
  cov_mat = zeros(n_sensor_locs);

  ## generate covariance matrix cov_mat

  for(i = 1:n_sensor_locs)
    for j=1:i
      dist_ij = dist([sensor_locs(1,i),sensor_locs(2,i)],[sensor_locs(1,j),sensor_locs(2,j)]);
      cov_mat(i,j) = covFun(dist_ij);
      cov_mat(j,i) = cov_mat(i,j);
     end
  end


  cov_krig = zeros(n_sensor_locs,1);
  est = zeros(n_kriging_locs,1);
  var = zeros(n_kriging_locs,1);

  ## estimate krigged mean and variance
  ## based on kriging maths
  
  for i = 1:n_kriging_locs
    kriging_loc = kriging_locs(1:2,i);
    for j = 1:n_sensor_locs
      dists = dist(kriging_loc',[sensor_locs(1,j),sensor_locs(2,j)]);
      cov_krig(j) = covFun(dists);
    end
    lambda = cov_mat^-1*cov_krig;
    est(i) = sensor_data*lambda;
    var(i) = cov_mat(1,1)-cov_krig'*lambda;
  end
end