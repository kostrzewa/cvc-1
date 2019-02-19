#define MAIN_PROGRAM
#include "global.h"
#undef MAIN_PROGRAM


#include "DependencyGraph.hpp"
#include "DependencyFulfilling.hpp"
#include "enums.hpp"
#include "types.h"
#include "Core.hpp"
#include "Logger.hpp"

#include <iostream>
#include <utility>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <memory>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/make_shared.hpp>
#include <fstream>
#include <boost/graph/graphml.hpp>

using namespace cvc;
using namespace boost;

int main(int argc, char ** argv)
{
  Core core(argc,argv,"boost_graph");
  Logger logger(0,0,std::cout);
  DepGraph g;

  std::vector<mom_t> in_momenta;
  std::vector<mom_t> out_momenta;

  for(int mom_x = -3; mom_x <= 3; mom_x++){
    for(int mom_y = -3; mom_y <= 3; mom_y++){
      for(int mom_z = -3; mom_z <= 3; mom_z++){
        if( mom_x*mom_x + mom_y*mom_y + mom_z*mom_z <= 1 ){
          in_momenta.push_back( mom_t{ mom_x, mom_y, mom_z } );
          out_momenta.push_back( mom_t{ mom_x, mom_y, mom_z } );
        }
      }
    }
  }

  int src_ts = 12;
  int dt = 12;

  for( auto const & in_mom : in_momenta ){
    double dummy_fwd_prop = in_mom.x;
    double dummy_bwd_prop = 0.0;

    char bwdpropname[200];
    snprintf(bwdpropname,
             200,
             "u/t23/g5/px0py0pz0");

    char fwdpropname[200];
    snprintf(fwdpropname,
             200,
             "u/t23/g%d/px%dpy%dpz%d",
             5,
             in_mom.x, in_mom.y, in_mom.z);

    for( auto const & out_mom : out_momenta ){
      mom_t seqmom = {-(out_mom.x+in_mom.x), 
                      -(out_mom.y+in_mom.y),
                      -(out_mom.z+in_mom.z)};

      char seqsrcname[200];
      snprintf(seqsrcname,
               200,
               "u/t%d/dt%d/gf%d/pfx%dpfy%dpfz%d",
               src_ts, dt,
               5, out_mom.x, out_mom.y, out_mom.z);

      Vertex seqsrcvertex = add_vertex(seqsrcname, g);
      g[seqsrcvertex].fulfill.reset( new SeqSourceFulfill(dt, out_mom, bwdpropname) );
      g[seqsrcvertex].independent = true;

      char seqpropname[200];
      snprintf(seqpropname,
               200,
               "sdu/gf%d/pfx%dpfy%dpfz%d/gi%d/pix%dpiy%dpiz%d",
               5, out_mom.x, out_mom.y, out_mom.z,
               5, in_mom.x, in_mom.y, in_mom.z);

      Vertex seqpropvertex = add_vertex(seqpropname, g);
      g[seqpropvertex].fulfill.reset( new PropFulfill("d", seqsrcname) );
      add_unique_edge(seqpropvertex, seqsrcvertex, g);

      char corrname[200];
      snprintf(corrname,
               200,
               "sdu+-g-u/gf%d/pfx%dpfy%dpfz%d/"
               "gc%d/"
               "gi%d/pix%dpiy%dpiz%d",
               5, out_mom.x, out_mom.y, out_mom.z,
               0,
               5, in_mom.x, in_mom.y, in_mom.z);
      
      Vertex corrvertex = add_vertex(corrname, g);
      g[corrvertex].fulfill.reset( new CorrFulfill(fwdpropname, seqpropname, seqmom, 0 ) ); 
      add_unique_edge(corrvertex, seqpropvertex, g);

      // for derivative operators we don't have any momentum transfer
      if( (in_mom.x == -out_mom.x && in_mom.y == -out_mom.y && in_mom.z == -out_mom.z) ||
          ( (in_mom.x == 0 && in_mom.x == out_mom.x) && 
            (in_mom.y == 0 && in_mom.y == out_mom.y) && 
            (in_mom.z == 0 && in_mom.z == out_mom.z) ) ){

        for( int dim1 : {DIM_T, DIM_X, DIM_Y, DIM_Z} ){
          for( int dir1 : {DIR_FWD, DIR_BWD}  ){
            for( int dim2 : {DIM_T, DIM_X, DIM_Y, DIM_Z} ){
              for( int dir2 : {DIR_FWD, DIR_BWD} ){
                for( int gc : { 0, 1, 2, 3, 4 } ){
                  char d1name[10];
                  snprintf(d1name,
                           10,
                           "d1_%c%c",
                           latDim_names[dim1],
                           shift_dir_names[dir1]);
                  char d2name[20];
                  snprintf(d2name,
                           20,
                           "d2_%c%c/%s",
                           latDim_names[dim2],
                           shift_dir_names[dir2],
                           d1name);

                  char Dpropname[200];
                  char DDpropname[200];
                  snprintf(Dpropname,
                           200,
                           "Du/%s/pix%dpiy%dpiz%d",
                           d1name,
                           in_mom.x, in_mom.y, in_mom.z);
                  snprintf(DDpropname,
                           200,
                           "DDu/%s/pix%dpiy%dpiz%d",
                           d2name,
                           in_mom.x, in_mom.y, in_mom.z);
        
                  Vertex Dpropvertex = add_vertex(Dpropname, g);
                  g[Dpropvertex].fulfill.reset( new CovDevFulfill( fwdpropname, dir1, dim1 ) );
                  g[Dpropvertex].independent = true; 

                  Vertex DDpropvertex = add_vertex(DDpropname, g);
                  add_unique_edge(DDpropvertex, Dpropvertex, g);
                  g[DDpropvertex].fulfill.reset( new CovDevFulfill( Dpropname, dir2, dim2 ) );
                  

                  char Dcorrname[200];
                  snprintf(Dcorrname,
                           200,
                           "sdu+-g-Du/gf%d/pfx%dpfy%dpfz%d/"
                           "gc%d/%s/"
                           "gi%d/pix%dpiy%dpiz%d",
                           5, out_mom.x, out_mom.y, out_mom.z,
                           gc, d1name,
                           5, in_mom.x, in_mom.y, in_mom.z);
                  
                  // even if we add this a million times, the vertex will be unique
                  Vertex Dcorrvertex = add_vertex(Dcorrname,g);
                  add_unique_edge(Dcorrvertex, Dpropvertex,g);
                  add_unique_edge(Dcorrvertex, seqpropvertex, g);
                  g[Dcorrvertex].fulfill.reset( new CorrFulfill(Dpropname, seqpropname, seqmom, gc ) ); 

                  char DDcorrname[200];
                  snprintf(DDcorrname,
                           200,
                           "sdu+-g-DDu/gf%d/pfx%dpfy%dpfz%d/"
                           "gc%d/"
                           "%s/"
                           "gi%d/pix%dpiy%dpiz%d",
                           5, out_mom.x, out_mom.y, out_mom.z,
                           gc, 
                           d2name,
                           5, in_mom.x, in_mom.y, in_mom.z);
                  Vertex DDcorrvertex = add_vertex(DDcorrname, g);
                  add_unique_edge(DDcorrvertex, DDpropvertex, g);
                  add_unique_edge(DDcorrvertex, seqpropvertex, g);
                  g[DDcorrvertex].fulfill.reset( new CorrFulfill(DDpropname, seqpropname, seqmom, gc ) );
                } // gc
              } // dir2
            } // dim2
          } // dir1
        } // dim1
      } // if(momenta)
    } // out_mom
  } // in_mom

  property_map<DepGraph, std::string VertexProperties::*>::type name_map = get(&VertexProperties::name, g);

  typedef graph_traits<DepGraph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  for(vp = vertices(g); vp.first != vp.second; ++vp.first){
    logger << name_map[*vp.first] << std::endl;
  }
  logger << std::endl;

  graph_traits<DepGraph>::edge_iterator ei, ei_end;
  graph_traits<DepGraph>::vertex_descriptor from, to;
  for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
    from = source(*ei, g);
    to = target(*ei, g);
    logger << "( " << g[from].name << "(" <<
      g[from].level << ")" << " -> " <<
      g[to].name << "(" << g[to].level << ")" <<
      " )" << std::endl;
  }
  logger << std::endl;

  boost::dynamic_properties dp;
  dp.property("name", name_map);
  std::fstream graphml_file("dependency_graph.graphml", std::ios_base::out);
  boost::write_graphml(graphml_file, g, dp, true);
  graphml_file.close();

  std::vector<int> component(num_vertices(g));
  int num = connected_components(g, &component[0]);
  std::vector<int>::size_type i;
  logger << "Total number of components: " << num << std::endl;
  for(i = 0; i != component.size(); ++i){
    g[i].component = component[i];
    logger << "Vertex " << name_map[i] << " is in component " << component[i];
    logger << " also stored " << g[i].component << std::endl;
  }
  logger << std::endl;

  std::vector<ComponentGraph> component_subgraphs(connected_components_subgraphs(g));
  for(size_t i_component = 0; i_component < component_subgraphs.size(); ++i_component){
    logger << std::endl << std::endl;
    logger << "Working on component " << i_component << std::endl;

    for( auto e : make_iterator_range(edges(component_subgraphs[i_component]))){
      std::cout << g[source(e,component_subgraphs[i_component])].name << " --> " << 
        g[target(e,component_subgraphs[i_component])].name << std::endl;
    }
    logger << std::endl;

    for( auto v : make_iterator_range(vertices(component_subgraphs[i_component]))){
      if( !g[v].fulfilled )
        descend_and_fulfill<DepGraph>( v, g );
    }
  }
  return 0;
}
