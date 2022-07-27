fitts_pointing_exgaussian_model = """

/*The data block for the pooled model.*/

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

/*Parameters, transformed parameters, and model blocks for the pooled model.*/

parameters {
    //Pooled Parameters
    real a;
    real<lower=0> b;
    real c;
    real<lower=0> d;
    real<lower=0> k;
}
transformed parameters {
    //Parameters used to make connections between mu, sigma and alpha, beta, lambda.
    real A = a;
    real B = b + 1 / k;
    real C = c;
    real D = d + square(1 / k);
    
    //The mean expressed by Fitts' law
    vector[number_of_conditions] mu = A + B * ID_list;
    //The standard deviation expressed by the quadratic variance model
    vector[number_of_conditions] sigma = sqrt(C + D * ID_square_list);
    
    //Parameter alpha of exGaussian model expressed by a, b, and ID
    vector[number_of_conditions] alpha = a + b * ID_list;
    //Parameter beta of exGaussian model expressed by c, d, and ID square
    vector[number_of_conditions] beta = sqrt(c + d * ID_square_list);
    //Parameter lambda of exGaussian model expressed by k and ID
    vector[number_of_conditions] lambda = k ./ ID_list;
}
model {
    a ~ normal(60.9, 100);
    b ~ normal(216.8, 100);
    c ~ normal(110.4, 10000);
    d ~ normal(716.7, 10000);
    k ~ normal(0.009, 100);
    for(i in 1 : number_of_conditions){
        movement_time[starts[i] + 1 : ends[i] + 1] ~ exp_mod_normal(alpha[i], beta[i], lambda[i]);
    }
}

/*The generated quantities block for the pooled model.*/
generated quantities {
    vector[number_of_data] log_likelihood;
    for(i in 1 : number_of_conditions){
        for(j in starts[i] + 1 : ends[i] + 1){
            log_likelihood[j]  = exp_mod_normal_lpdf(movement_time[j] | alpha[i], beta[i], lambda[i]);
        }
    }
}
"""