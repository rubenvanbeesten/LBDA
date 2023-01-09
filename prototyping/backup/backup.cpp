/*
Prototype for dealing with risk measures in LBDA

Inputs:

lambdas:            vector of weights for each cvar
betas:              vector of cvar parameters (in increasing order) #TODO: make check of this fact
scenario_values:    vector of values corresponding to each value

Outputs: 

weights:            vector of weights per scenario
phi:                risk spectrum function
sorted_idx:         sorted indices corresponding to scenario_values
rho:                risk measure applied to scenario_values
*/



// Includes
#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

// Note: no using namespace std! (also not in original LBDA code)


// 1. SORTING

template <typename T> // this indicates that "T" is a typename

// Function that sorts and stores indices 
std::vector<std::size_t> sort_indices(const std::vector<T> &v) {

    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    std::stable_sort(idx.begin(), idx.end(),
        [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});

    return idx;
}

// Function that tests sorting
void test_sorting() {
    std::cout << "Testing sorting." << std::endl;

    // input vector
    std::vector<double> v = {1.0, 3.0, 4.0, 6.0, 2.0, 5.0};

    std::cout << "Original v: " << std::endl;
    for (int i = 0; i < v.size(); i++){
        std::cout << i << ": " << v.at(i) << std::endl;
    }

    std::cout << "Sort indices..." << std::endl;

    std::vector<int> sorted_idx(v.size());
    int idx_counter = 0;
    for (auto i: sort_indices(v)) { 
        sorted_idx.at(idx_counter) = i;
        idx_counter++;
    }

    std::cout << "Output sorted indices and sorted vector:" << std::endl;

    for (int i=0; i < sorted_idx.size(); i++){
        std::cout << sorted_idx.at(i) << ": " << v.at(sorted_idx.at(i)) << std::endl;
    }
}


// 2. RISK MEASURE

// 2a. Check risk measure

// Perform sensibility checks on risk measure
void check_risk_measure(std::vector<double> lambdas, std::vector<double> betas) {
    // Throw error if:
    // - lambdas and betas are not same length
    // - betas not increasing
    // - lambdas don't add to one
    // - betas not in [0,1)
    if(lambdas.size() != betas.size()){
        throw std::runtime_error("lambdas and betas should have same length");
    } 
    double sum_lambdas = 0.0;
    for(int i=0; i < lambdas.size(); i++){
        sum_lambdas += lambdas.at(i);
        if(betas.at(i) < 0.0 || betas.at(i) >= 1.0){
            throw std::runtime_error("betas should be in interval [0, 1)");
        }
        if(i > 0){
            if(betas.at(i) <= betas.at(i-1)){
                throw std::runtime_error("betas should be strictly increasing");
            }
        }
    }
    if (std::abs(sum_lambdas - 1.0) > 0.001) {
        throw std::runtime_error("lambdas should add to one");
    }
}


// 2b. Phi 

// Compute phi (risk spectrum)
// as a function of p in (0,1) and the risk measure
double phi(double p, std::vector<double> lambdas, std::vector<double> betas)
{
    // Sanity checks on the risk measure
    check_risk_measure(lambdas, betas);
    
    // Compute phi
    double res = 0.0; // initialize at zero
    for(int i=0; i < lambdas.size(); i++){
        if(p >= betas[i]){
            res += lambdas.at(i) / (1.0 - betas.at(i)); //iteratively add pieces
        }
    }

    return res;
}

// Test phi function
void test_phi() {
    std::cout << "Testing phi." << std::endl;

    std::vector<double> lambdas = {0.25, 0.25, 0.25, 0.25};
    std::vector<double> betas = {0.0, 0.5, 0.9, 0.98};

    for (double p = 0.0; p < 1.0; p += 0.01) {
        double cur_phi = phi(p, lambdas, betas);
        std::cout << "phi(" << p << ") = " << cur_phi << std::endl;
    }

    std::cout << "Finished testing phi." << std::endl;

}


// 2c. Compute probability weights

// Function that computes vector of weights based on risk measure
// weights are increasing
std::vector<double> compute_weights(int num_scenarios, std::vector<double> lambdas, std::vector<double> betas)
{
    // Sensibility check on risk measure
    check_risk_measure(lambdas, betas);

    // initialize vector of weights
    std::vector<double> weights(num_scenarios, 0.0);

    // compute weights
    for(int s=0; s<num_scenarios; s++){

        double left = (1.0*s)/num_scenarios; //left endpoint of interval
        double right = left + 1.0/num_scenarios; //right endpoint of interval

        //initialize with starting height times width
        double cur_weight = phi(left, lambdas, betas) * (right - left);

        //next, add any additional blocks
        for(int i = 0; i < lambdas.size(); i++){
            if((betas.at(i) > left) && (betas.at(i) < right)){                                  //if betas[i] == left, then already included in phi(left)
                cur_weight += (right - betas.at(i)) * (lambdas.at(i) / (1.0 - betas.at(i)));    //add piece of this step that falls within the interval
            }
        }

        weights.at(s) = cur_weight;
    }

    return weights;

}

// Test compute_weights

void test_weights() {
    std::cout << "Testing compute_weights" << std::endl;

    // risk measure
    std::vector<double> lambdas = {0.25, 0.25, 0.25, 0.25};
    std::vector<double> betas = {0.0, 0.5, 0.9, 0.98};

    // scenarios
    int num_scenarios = 100;

    // compute weights
    std::vector<double> weights = compute_weights(num_scenarios, lambdas, betas);

    // print weights
    std::cout << "Printing weights:" << std::endl;
    for(int i=0; i < weights.size(); i++){
        std::cout << "weights[" << i << "] = " << weights.at(i) << std::endl;
    }

    std::cout << "Done testing compute_weights" << std::endl;

}

// 2d. Compute rho

double compute_rho(std::vector<double> scenario_values, std::vector<double> lambdas, std::vector<double> betas){
    // extract number of scenarios
    int num_scenarios = scenario_values.size();

    // compute weights
    std::vector<double> weights = compute_weights(num_scenarios, lambdas, betas);

    // sort scenario values
    std::vector<int> sorted_idx(scenario_values.size()); // this will contain the sorted indices
    int idx_counter = 0;
    for (auto i: sort_indices(scenario_values)) { 
        sorted_idx.at(idx_counter) = i;
        idx_counter++;
    }

    // compute rho
    double res = 0.0; // initialize result
    for(int s=0; s < scenario_values.size(); s++){
        res += weights.at(s) * scenario_values.at(sorted_idx.at(s));
    }

    return res;
}

// Testing rho
void test_rho() {
    std::cout << "Testing compute_rho..." << std::endl;

    // define risk measure
    std::vector<double> lambdas = {0.25, 0.25, 0.25, 0.25};
    std::vector<double> betas = {0.0, 0.5, 0.9, 0.98};

    // provide scenario values
    std::vector<double> scenario_values = {4.0, 5.0, 2.0, 3.0, 1.0, 7.0, 8.0, 9.0, 0.0, 6.0};

    // compute rho
    double rho = compute_rho(scenario_values, lambdas, betas);

    // print rho
    std::cout << "rho(scenario_values) = " << rho << std::endl;

    std::cout << std::endl << "Done with prototype example." << std::endl;
}



// 3. Testing everything

void prototype_example() {
    // 1. INPUT

    // define risk measure
    std::vector<double> lambdas = {0.25, 0.25, 0.25, 0.25};
    std::vector<double> betas = {0.0, 0.5, 0.9, 0.98};

    // provide scenario values
    std::vector<double> scenario_values = {4.0, 5.0, 2.0, 3.0, 1.0, 7.0, 8.0, 9.0, 0.0, 6.0};
    int num_scenarios = scenario_values.size();

    // 2. COMPUTE WEIGHTS

    // compute weights
    std::vector<double> weights = compute_weights(num_scenarios, lambdas, betas);

    // 3. SORT SCENARIO VALUES

    // sort scenario values
    std::vector<int> sorted_idx(scenario_values.size()); // this will contain the sorted indices
    int idx_counter = 0;
    for (auto i: sort_indices(scenario_values)) { 
        sorted_idx.at(idx_counter) = i;
        idx_counter++;
    }

    // 4. COMPUTE RHO

    // compute rho(scen_values)
    double rho = compute_rho(scenario_values, lambdas, betas);


    // 5. OUTPUT

    // print weights
    std::cout << "Weights:" << std::endl;
    for(int i=0; i < weights.size(); i++){
        std::cout << "weights[" << i << "] = " << weights.at(i) << std::endl;
    }
    std::cout << std::endl;

    // print sorted values
    std::cout << "Sorted scenario_values:" << std::endl;
    for (int i=0; i < sorted_idx.size(); i++){
        std::cout << sorted_idx.at(i) << ": " << scenario_values.at(sorted_idx.at(i)) << std::endl;
    }
    std::cout << std::endl;

    // print rho
    std::cout << "rho(scenario_values) = " << rho << std::endl;

    std::cout << std::endl << "Done with prototype example." << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////////////




// 4. Add everything together in a class

class RiskMeasure
{
protected:
    std::vector<double> lambdas;
    std::vector<double> betas;

public:
    // constructor
    RiskMeasure(std::vector<double> lambdas_new, std::vector<double> betas_new);

    // destructor
    ~RiskMeasure();

    // compute phi
    double phi(double p);

    // compute vector of weights
    std::vector<double> compute_weights(int num_scenarios);

    // compute rho
    double compute_rho(std::vector<double> scenario_values);

};

// constructor
RiskMeasure::RiskMeasure(std::vector<double> lambdas_new, std::vector<double> betas_new) {
    
    // Set values
    lambdas = lambdas_new;
    betas = betas_new;

    // Sensibility checks: throw error if
    // - lambdas and betas are not same length
    // - betas not increasing
    // - lambdas don't add to one
    // - betas not in [0,1)

    if(lambdas.size() != betas.size()){
        throw std::runtime_error("lambdas and betas should have same length");
    } 
    double sum_lambdas = 0.0;
    for(int i=0; i < lambdas.size(); i++){
        sum_lambdas += lambdas.at(i);
        if(betas.at(i) < 0.0 || betas.at(i) >= 1.0){
            throw std::runtime_error("betas should be in interval [0, 1)");
        }
        if(i > 0){
            if(betas.at(i) <= betas.at(i-1)){
                throw std::runtime_error("betas should be strictly increasing");
            }
        }
    }
    if (std::abs(sum_lambdas - 1.0) > 0.001) {
        throw std::runtime_error("lambdas should add to one");
    }
}

// destructor
RiskMeasure::~RiskMeasure(void) {
    // do nothing
}


// Compute phi (risk spectrum)
// as a function of p in (0,1) and the risk measure
double RiskMeasure::phi(double p)
{  
    double res = 0.0; // initialize at zero
    for(int i=0; i < lambdas.size(); i++){
        if(p >= betas[i]){
            res += lambdas.at(i) / (1.0 - betas.at(i)); //iteratively add pieces
        }
    }

    return res;
}

// Compute vector of weights (weights are increasing)
std::vector<double> RiskMeasure::compute_weights(int num_scenarios)
{
    // initialize vector of weights
    std::vector<double> weights(num_scenarios, 0.0);

    // compute weights
    for(int s=0; s<num_scenarios; s++){

        double left = (1.0*s)/num_scenarios; //left endpoint of interval
        double right = left + 1.0/num_scenarios; //right endpoint of interval

        //initialize with starting height times width
        double cur_weight = RiskMeasure::phi(left) * (right - left);

        //next, add any additional blocks
        for(int i = 0; i < lambdas.size(); i++){
            if((betas.at(i) > left) && (betas.at(i) < right)){                                  //if betas[i] == left, then already included in phi(left)
                cur_weight += (right - betas.at(i)) * (lambdas.at(i) / (1.0 - betas.at(i)));    //add piece of this step that falls within the interval
            }
        }

        weights.at(s) = cur_weight;
    }

    return weights;
}

// 2d. Compute rho

double RiskMeasure::compute_rho(std::vector<double> scenario_values){
    // extract number of scenarios
    int num_scenarios = scenario_values.size();

    // compute weights
    std::vector<double> weights = RiskMeasure::compute_weights(num_scenarios);

    // sort scenario values
    std::vector<int> sorted_idx(scenario_values.size()); // this will contain the sorted indices
    int idx_counter = 0;
    for (auto i: sort_indices(scenario_values)) { 
        sorted_idx.at(idx_counter) = i;
        idx_counter++;
    }

    // compute rho
    double res = 0.0; // initialize result
    for(int s=0; s < scenario_values.size(); s++){
        res += weights.at(s) * scenario_values.at(sorted_idx.at(s));
    }

    return res;
}



// Main function
int main() {
    try{
        //test_sorting(); 
        //test_phi();  
        //test_weights();
        //test_rho();
        //prototype_example();

        // construct risk measure
        std::vector<double> lambdas = {0.5, 0.5};
        std::vector<double> betas = {0.0, 0.9};
        RiskMeasure risk_measure(lambdas, betas);

        // compute phi
        std::cout << "phi(0.5) = " << risk_measure.phi(0.5) << std::endl;
        std::cout << "phi(0.95) = " << risk_measure.phi(0.95) << std::endl;
        std::cout << std::endl;    

        // compute vector of weights
        int num_scenarios = 10;
        std::vector<double> weights = risk_measure.compute_weights(num_scenarios);
        std::cout << "Vector of weights:" << std::endl;
        for(int i = 0; i < num_scenarios; i++){
            std::cout << "w[" << i << "] = " << weights[i] << std::endl;
        }
        std::cout << std::endl;          

        // compute rho
        std::vector<double> scenario_values = {0.0, 9.0, 8.0, 7.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        double rho = risk_measure.compute_rho(scenario_values);
        std::cout << "rho = " << rho << std::endl;
        
        std::cout << std::endl << "Finished testing." << std::endl;


        return 0;
    }



    catch(std::exception& e){

        std::cerr << "Error: " << e.what() << '\n';

        // do nothing
        return -1;
    }
}