function[ cred, pbayesmean, alpha_post, beta_post] = beta_binomial_cred_interval(alpha, beta, correct, p)

% this function computes the p*100% credible interval from the bernoulli
% vector correct using beta hyperparmeter values of alpha and beta

% uniform prior 
% alpha = 1; 
% beta = 1; 

% compute posterior hyperparameters
alpha_post = alpha + sum(correct); 
beta_post = beta + length(correct) - sum(correct);

% compute the posterior mean on success probability 
pbayesmean = alpha_post/(alpha_post + beta_post); 

% compute the credible interval 
plow = (1 -p)/2; % lower probability 
pup = p + plow; % upper probability 
lower = betainv(plow,alpha_post,beta_post);
upper = betainv(pup,alpha_post,beta_post);

cred = [lower, upper]; 


% %% Solve for optimal alpha, beta given chance level 
% 
% c = 1/174; 
% beta = 3; 
% alpha = (beta*c - 2*c + 1)/(1-c); % mode at c 
% 
% alpha;
% beta;
