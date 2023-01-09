// This file contains the risk measure class, plus a sorting function

#include "risk_measure.h"

/////////////////////
// RISK MEASURE CLASS
/////////////////////

// Constructor from vectors
RiskMeasure::RiskMeasure(std::vector<double> lambdas_new, std::vector<double> betas_new) {
    
    // Set values
    lambdas = lambdas_new;
    betas = betas_new;

    // do a sensibility check
    RiskMeasure::sensibilityCheck();
    
}

// Destructor
RiskMeasure::~RiskMeasure(void) {
    // do nothing
}

// sensibility check
void RiskMeasure::sensibilityCheck(void)
{
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

// Compute rho

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

/////////////////////////////////////////////////////////

/////////////////////
// SORTING FUNCTION
/////////////////////


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

//////////////////////////////
// STRING TO DOUBLE FUNCTION
//////////////////////////////


// Function that converts string to vector of double
// - input of the form: 0.2,0.3,0.5 (no spaces!)
std::vector<double> str_to_double_vec(std::string s){
    // create stringstream from s
    std::stringstream ss(s);

    // create vector to be returned
    std::vector<double> result;

    // create temporary string
    std::string tmp;

    // iteratively push to vector "words"
    while(getline(ss, tmp, ',')){
        result.push_back(std::stod(tmp));
    }

    return result;

}
