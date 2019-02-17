#include <yaml-cpp/yaml.h>
#include <iostream>

#include "meta_types.hpp"
#include "yaml_parsers.hpp"

namespace cvc {
namespace yaml {

void enter_node(const YAML::Node &node, 
                const unsigned int depth,
                MetaCollection & metas,
                const bool verbose){
  YAML::NodeType::value type = node.Type();
  std::string indent( 2*(size_t)depth, ' ');
  switch(type){
    case YAML::NodeType::Scalar:
      if(verbose){
        std::cout << node;
      }
      break;
    case YAML::NodeType::Sequence:
      for(unsigned int i = 0; i < node.size(); ++i){
        if(verbose && depth <= 2){
          std::cout << "[ ";
        }
        const YAML::Node & subnode = node[i];
        enter_node(subnode, depth+1, metas, verbose);
        if(verbose){
          if( depth <= 2 ){
            std::cout << " ]";
          } else if( i < node.size()-1 ) {
            std::cout << ",";
          }
        }
      }
      break;
    case YAML::NodeType::Map:
      for(YAML::const_iterator it = node.begin(); it != node.end(); ++it){
        if(verbose){
          if(depth != 0){ std::cout << std::endl << indent; }
          std::cout << it->first << ": ";
        }
        if( it->first.as<std::string>() == "MomentumList" ){
          construct_momentum_list(it->second, verbose, metas.mom_lists );
        } else if ( it->first.as<std::string>() == "TimeSlicePropagator" ){
          construct_time_slice_propagator(it->second, verbose, metas.mom_lists, metas.ts_props );
        } else if ( it->first.as<std::string>() == "OetMesonTwoPointFunction" ){
          construct_oet_meson_two_point_function(it->second, verbose, metas.mom_lists, metas.ts_props, metas.g); 
        } else {
          char msg[200];
          snprintf(msg, 200,
                   "%s is not a valid Object name\n",
                   it->first.as<std::string>().c_str());
          throw( std::invalid_argument(msg) );
        }
      }
    case YAML::NodeType::Null:
      if(verbose)
        std::cout << std::endl;
      break;
    default:
      break;
  }
  for(const auto & prop : metas.ts_props ){
    std::cout << prop.first << std::endl;
  }
}

} // namespace(yaml)
} // namespace(cvc)