%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%' Simulate paths of GBM
%'
%' Simulate from p(X_t|X_{t-1})
%' Use log-normal transition density specified by the model
%' @param x0 is the starting values
%' @param dt is the step size
%' @param model contains all the other parameters
%' @export

function newX = sim_gbm( x0, model, dt )

if (nargin == 2)
    dt = model.dt;
end

len = size(x0,1);

newX = x0;

for j = 1:model.dim
   newX(:,j) = x0(:,j).*exp( normrnd(0,1,[len,1])*model.sigma(j)*sqrt(dt) + (model.r- model.div- model.sigma(j)^2/2)*dt );
end

end