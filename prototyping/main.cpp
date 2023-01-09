#include "risk_measure.h"

// Main function
int main() {
    try{
       
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