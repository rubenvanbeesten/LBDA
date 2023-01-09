#include <iostream>
#include <vector>
#include <sstream>

// function that translates a string to a vector of doubles
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


int main(){

    std::string s = "0.1,0.2,0.3,0.4";

    auto vec = str_to_double_vec(s);

    for(int i=0; i<vec.size(); i++){
        std::cout << vec.at(i) << std::endl;
    }
}
