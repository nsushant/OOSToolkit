#include "data_access_lib.hpp"

arma::uvec find_idxs_of_match(std::vector<std::string> s, std::string value_to_match){

    std::string target = value_to_match;
    std::vector<int> idxs;

    for (int i = 0; i < (int)s.size() ; ++i) {

        if (s[i] == target) {
            idxs.push_back(i);
        }
    }

    // Convert to arma::uvec
    arma::uvec arma_indices(idxs.size());

    for (size_t i = 0; i < idxs.size(); ++i) {
        arma_indices(i) = static_cast<arma::uword>(idxs[i]);
    }


    return arma_indices;



    }
