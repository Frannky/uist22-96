fitts_pointing_exgaussian_model = """

/*The data block for the hierarchical model.*/

data {
    //Number of individuals
    int number_of_individuals;
    //Number of conditions (user, task) pair
    int number_of_conditions;
    //Number of data points
    int number_of_data;
    //Index of difficulty (ID) list
    vector[number_of_conditions] ID_list;
    //Index of difficulty squared (ID^2) list
    vector[number_of_conditions] ID_square_list;
    //individual index of each condition
    int<lower=1, upper=number_of_individuals> individual_index[number_of_conditions];
    //Movement time list for all data points
    vector[number_of_data] movement_time;
    //Starting index of each condition in movement time list
    int starts[number_of_conditions];
    //Ending index of each condition in movement time list
    int ends[number_of_conditions];
}

/*Parameters, transformed parameters, and model blocks for the hierarchical model.*/

parameters {
    //The mean of the hyperprior of parameter a, b, c, d, k
    real mu_a;
    real<lower=0> mu_b;
    real mu_c;
    real<lower=0> mu_d;
    real<lower=0> mu_k;
    //The standard deviation of the hyperprior of parameter a, b, c, d, k
    real<lower=0> tau_a;
    real<lower=0> tau_b;
    real<lower=0> tau_c;
    real<lower=0> tau_d;
    real<lower=0> tau_k;
    //eta_a, eta_b, eta_c, eta_d, and eta_k are standard normal distributions
    vector[number_of_individuals] eta_a;
    vector<lower=0>[number_of_individuals] eta_b;
    vector<lower=0>[number_of_individuals] eta_c;
    vector<lower=0>[number_of_individuals] eta_d;
    vector<lower=0>[number_of_individuals] eta_k;
}
transformed parameters {
    //The mean (mu) and standard deviation (sigma) of the distribution
    vector[number_of_conditions] mu;
    vector[number_of_conditions] sigma;
    
    //The location parameter alpha, scale parameter beta, and shape parameter lambda of the exGaussian distribution
    vector[number_of_conditions] alpha;
    vector[number_of_conditions] beta;
    vector[number_of_conditions] lambda;

    //Parameters associated with each user
    vector[number_of_individuals] a = mu_a + tau_a * eta_a;
    vector[number_of_individuals] b = mu_b + tau_b * eta_b;
    vector[number_of_individuals] c = mu_c + tau_c * eta_c;
    vector[number_of_individuals] d = mu_d + tau_d * eta_d;
    vector[number_of_individuals] k = mu_k + tau_k * eta_k;
    
    //Parameters used to make connections between mu, sigma and alpha, beta, lambda.
    vector[number_of_individuals] A = a;
    vector[number_of_individuals] B = b + 1 ./ k;
    vector[number_of_individuals] C = c;
    vector[number_of_individuals] D = d + 1 ./ square(k);

    for (i in 1:number_of_conditions){
        //The mean expressed by Fitts' law
        mu[i] = A[individual_index[i]] + B[individual_index[i]] * ID_list[i];
        //The standard deviation expressed by the quadratic variance model
        sigma[i] = sqrt(C[individual_index[i]] + D[individual_index[i]] * ID_square_list[i]);
    }

    for (i in 1:number_of_conditions) {
        //Parameter alpha of exGaussian model expressed by a, b, and ID
        alpha[i] = b[individual_index[i]] * ID_list[i] + a[individual_index[i]];
        //Parameter beta of exGaussian model expressed by c, d, and ID square
        beta[i] = sqrt(d[individual_index[i]] * ID_square_list[i] + c[individual_index[i]]);
        //Parameter lambda of exGaussian model expressed by k and ID
        lambda[i] = k[individual_index[i]] / ID_list[i];
    }
}
model {
    eta_a ~ normal(0, 1);
    eta_b ~ normal(0, 1);
    eta_c ~ normal(0, 1);
    eta_d ~ normal(0, 1);
    eta_k ~ normal(0, 1);
    mu_a ~ normal(18.2, 100);
    tau_a ~ cauchy(32.1, 100);
    mu_b ~ normal(20.2, 100);
    tau_b ~ cauchy(59.4, 100);
    mu_c ~ normal(1.8, 10000);
    tau_c ~ cauchy(3.9, 10000);
    mu_d ~ normal(19.8, 10000);
    tau_d ~ cauchy(48.6, 10000);
    mu_k ~ normal(0.02, 100);
    tau_k ~ cauchy(0.01, 100);
    for(i in 1 : number_of_conditions){
        movement_time[starts[i] + 1 : ends[i] + 1] ~ exp_mod_normal(alpha[i], beta[i], lambda[i]);
    }
}

/*The generated quantities block for the hierarchical model.*/
generated quantities {
    vector[number_of_data] log_likelihood;
    for(i in 1 : number_of_conditions){
        for(j in starts[i] + 1 : ends[i] + 1){
            log_likelihood[j]  = exp_mod_normal_lpdf(movement_time[j] | alpha[i], beta[i], lambda[i]);
        }
    }
}
"""