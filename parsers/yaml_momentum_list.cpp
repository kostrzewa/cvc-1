#include <yaml-cpp/yaml.h>
#include <vector>
#include <map>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cmath>

#include "algorithms.hpp"
#include "types.h"

namespace cvc {
namespace yaml {

std::vector< mom_t > parse_momentum_list(const YAML::Node & node)
{
  if( node.Type() != YAML::NodeType::Sequence ){
    throw( std::invalid_argument("in parse_momentum_list, 'node' must be of type YAML::NodeType::Sequence\n") );
  }
  std::vector< mom_t > momenta;
  for(size_t i = 0; i < node.size(); ++i){
    mom_t mom{ node[i][0].as<int>(), node[i][1].as<int>(), node[i][2].as<int>() };
    momenta.push_back(mom);
  }
  momentum_compare_t comp;
  std::sort(momenta.begin(), momenta.end(), comp);
  return momenta;
}

std::vector< mom_t > psq_to_momentum_list(const YAML::Node & node)
{
  if( node.Type() != YAML::NodeType::Scalar ){
    throw( std::invalid_argument("in psq_to_momentum_list, 'node' must be of type YAML::NodeType::Scalar\n") );
  }
  const int psqmax = node.as<int>();
  const int pmax = static_cast<int>(sqrt( node.as<double>() )); 
  std::vector< mom_t > momenta;
  for( int px = -pmax; px <= pmax; ++px ){
    for( int py = -pmax; py <= pmax; ++py ){
      for( int pz = -pmax; pz <= pmax; ++pz ){
        if( (px*px + py*py + pz*pz) <= psqmax ){
          mom_t mom{ px, py, pz };
          momenta.push_back(mom);
        }
      }
    }
  }
  momentum_compare_t comp;
  std::sort(momenta.begin(), momenta.end(), comp);
  return momenta; 
}

void construct_momentum_list(const YAML::Node & node,
                             const bool verbose,
                             mom_lists_t & mom_lists )
{
  if( node.Type() != YAML::NodeType::Map ){
    throw( std::invalid_argument("in construct_momentum_list, 'node' must be of type YAML::NodeType::Map\n") );
  }
  if( !node["id"] || !(node["Psqmax"] || node["Plist"]) ){
    throw( std::invalid_argument("For 'MomentumList', the 'id' property and one of 'Psqmax' or 'Plist' must be defined!\n") );
  }

  std::string id;
  std::vector<mom_t> momentum_list;
  for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
    if(verbose){
      std::cout << "\n  " << it->first << ": " << it->second;
    }
    if( it->first.as<std::string>() == "id" ){
      id = it->second.as<std::string>();
    } else if( it->first.as<std::string>() == "Plist" ){
      momentum_list = parse_momentum_list( it->second );
    } else if( it->first.as<std::string>() == "Psqmax" ){
      momentum_list = psq_to_momentum_list( it->second );
    } else {
      char msg[200];
      snprintf(msg, 200,
               "%s is not a valid property for a MomentumList!\n",
               it->first.as<std::string>().c_str());
      throw( std::invalid_argument(msg) );
    }
  }
  if( verbose ){
    std::cout << std::endl << "   ";
    for( const auto & mom : momentum_list ){
      std::vector<int> mvec{mom.x, mom.y, mom.z};
      std::cout << "[";
      for(size_t i = 0; i < mvec.size(); ++i){
        std::cout << mvec[i];
        if( i < mvec.size()-1 ){
          std::cout << ",";
        } else {
          std::cout << "] ";
        }
      }
    }
    std::cout << std::endl;
  }
  mom_lists[id] = momentum_list;
}

} //namespace(yaml)
} //namespace(cvc)
