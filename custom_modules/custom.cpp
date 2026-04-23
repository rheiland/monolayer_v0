/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <algorithm>    // for std::remove

#include "./custom.h"

double time_90pct;  // time at which the outer cells reach 90% max width
double xpos_90pct;
bool cells11_flag = true;

double rest_length_factor = 1.0;

// double tmin=4;
// double tmax= 5;
double tmin=0;
double tmax= 1;

void create_cell_types( void )
{
    if (parameters.doubles.find_index("rest_length_factor") != -1)
	{
        rest_length_factor = parameters.doubles("rest_length_factor");
	}
    else
    {
        std::cout << "Error: rest_length_factor needs to be defined in user params\n" << std::endl;
        std::exit(-1);
    }
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
    cell_defaults.functions.cell_division_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	// cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_cell_rule; 
    cell_defaults.functions.update_velocity = custom_update_cell_velocity;
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    xpos_90pct = (0.9 * cell_defaults.phenotype.geometry.radius * 2 * 10.0) / 2.0;
    std::cout << "-------- xpos_90pct= " << xpos_90pct << std::endl;
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
    // First, decide if we're doing the 11 or 21 cell relaxation model
    cells11_flag = parameters.bools("cells11");

	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	

	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }


// int Cell::get_current_mechanics_voxel_index()
// {
// 	return current_mechanics_voxel_index;
// }


// void Cell::update_voxel_in_container()
// {
// 	// call the method from BioFVM_basic_agent to update microenvironment's voxel index
// 	update_voxel_index();
// 	// int temp_current_voxel_index;
// 	// Check to see if we need to remove agents that are pushed out of boundary
// 	// if(!get_container()->underlying_mesh.is_position_valid(position[0],position[1],position[2]))	
		
// 	if(updated_current_mechanics_voxel_index==-1)// updated_current_mechanics_voxel_index is updated in update_position
// 	{
// 		// check if this agent has a valid voxel index, if so, remove it from previous voxel
// 		if( get_current_mechanics_voxel_index() >= 0)
// 		{
// 			{get_container()->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());}
// 		}
// 		{get_container()->add_agent_to_outer_voxel(this);}
// 		// std::cout<<"cell out of boundary..."<< __LINE__<<" "<<ID<<std::endl;
// 		current_mechanics_voxel_index=-1;
// 		is_out_of_domain=true;
// 		is_active=false;
// 		return;
// 	}
	
// 	// temp_current_voxel_index= get_current_mechanics_voxel_index();
// 	// updated_current_mechanics_voxel_index=get_container()->underlying_mesh.nearest_voxel_index( position );
	
// 	// update mesh indices (if needed)
// 	if(updated_current_mechanics_voxel_index!= get_current_mechanics_voxel_index())
// 	{
// 		{
// 			container->remove_agent_from_voxel(this, get_current_mechanics_voxel_index());
// 			container->add_agent_to_voxel(this, updated_current_mechanics_voxel_index);
// 		}
// 		current_mechanics_voxel_index=updated_current_mechanics_voxel_index;
// 	}
	
// 	return; 
// }

void get_cells_in_mech_voxels()
{
    // iterate over all mech voxels

}

// This provides a "correction" for the pCell->state.neighbors determined in the core code.
// Specifically, it enforces a more constrained, distance measure of what is a neighbor cell.
// void generate_neighbor_list( Cell* pCell)
// {
//         // double tmin=479;
//         // double tmax=480;
//     // for debug prints
//     // static double tmin= 4.0;
//     // static double tmax= 1331;
//     // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//     //     std::cout << "-------- " <<__FUNCTION__ << ":  t= "<<PhysiCell_globals.current_time << ": pCell->ID= " << pCell->ID<< std::endl;

//     if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//     {
//         std::cout << "\n" << __FUNCTION__ << "  t= "<< PhysiCell_globals.current_time << ":  ~~~~~~~~~~~~~~~~~~~~   pCell->ID = " << pCell->ID << std::endl;
//     }

// 	pCell->state.neighbors.clear();
//     std::vector<PhysiCell::Cell *> true_nbrs;

// 	// ---------- 1) First check the neighbors in my current voxel --------------
//     if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//     {
//         // std::cout << "   my voxel center= " << pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center << std::endl;
//         std::cout << "   my current voxel index= " << pCell->get_current_mechanics_voxel_index() << std::endl;
//     }

// 	std::vector<Cell*>::iterator neighbor;
// 	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();

// 	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
// 	{
// 		// add_potentials_monolayer(pCell, *neighbor);

//         Cell* pN = *neighbor;
//         if (pCell->ID == pN->ID) continue;

//         if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//         {
//             std::cout << __FUNCTION__ << "--- (current voxel): ID=" << pCell->ID << ", pN->ID= " << pN->ID << std::endl;
//         }


//         bool nbr_exists = false;
//         for (const auto& nbr : pCell->state.neighbors)
//         {
//             if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//             {
//                 std::cout << __FUNCTION__ << "---------nbr->ID = " << nbr->ID << std::endl;
//             }
//             if (nbr->ID == pN->ID)
//             {
//                 if (PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//                     std::cout << __FUNCTION__ << "------------- same, so exit " << std::endl;
//                 nbr_exists = true;
//                 break;
//             }
//         }
//             if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//             {
//                 std::cout << __FUNCTION__ << "------ nbr_exists=  " << nbr_exists << std::endl;
//             }
//         if (!nbr_exists)
//         {
//             pCell->state.neighbors.push_back(*neighbor);
//             if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//             {
//                 std::cout << __FUNCTION__ << "------------- DNE, so push_back " << std::endl;
//                 for (const auto& nbr : pCell->state.neighbors)
//                     std::cout << __FUNCTION__ << "---------------- ID= " << pCell->ID << " and new nbr->ID = " << nbr->ID << std::endl;
//             }
//         }

//         if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//             std::cout << "   (current voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << std::endl;

//         // std::cout << __FUNCTION__ << "  t= "<<PhysiCell_globals.current_time << ": (current voxel) ID= " << pCell->ID<< ",  neighbor ID=" << pN->ID << ", pCell num_nbrs= "<< pCell->custom_data["num_nbrs"] << std::endl;
// 	}


// 	// ---------- 2) Second, check the cells in the surrounding voxels (mechanics voxel size) --------------
//     if (PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//         std::cout << "    2) ~~~~~~~~~~~~~   surrounding voxels" << std::endl;
// 	std::vector<int>::iterator neighbor_voxel_index;
// 	std::vector<int>::iterator neighbor_voxel_idx_end = 
// 		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();


// 	for( neighbor_voxel_index = 
// 		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
// 		neighbor_voxel_index != neighbor_voxel_idx_end; 
// 		++neighbor_voxel_index )
// 	{
//         if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//             std::cout << "   nbr voxel center= " << pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center << std::endl;

// 		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
// 			continue;

// 		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();

// 		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
// 		{
//             Cell* pN = *neighbor;
//             if (pCell->ID == pN->ID) continue;

//             bool nbr_exists = false;
//             for (const auto& nbr : pCell->state.neighbors)
//             {
//                 if (nbr->ID == pN->ID)
//                 {
//                     nbr_exists = true;

//                     if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//                     {
//                         std::cout << "       -- for ID=1, nbr_exists ID=" << nbr->ID << std::endl;
//                     }
//                     break;
//                 }
//             }
//             if (!nbr_exists)
//             {
//                 pCell->state.neighbors.push_back(*neighbor);

//                 if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//                 {
//                     std::cout << "       -- for ID=1, push_back ID=" << (*neighbor)->ID << std::endl;
//                 }
//             }
// 		}
// 	}

//     // if (PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//     // {
//     //     std::cout << " -- nbrs, minus duplicates:" << std::endl;
//     //     for (const auto& nbr : pCell->state.neighbors)
//     //         std::cout <<  nbr->ID << std::endl;
//     // }

//     // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//     // {
//     //     std::cout << " -------- ID 1 has prelim nbrs: " << std::endl;
//     //     for (const auto& nbr : pCell->state.neighbors)
//     //         std::cout <<  nbr->ID << std::endl;
//     // }


//     // -----------  3) Lastly, cull the nbrs to match a distance constraint  --------------
//     // if (PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//     if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//     {
//         std::cout << "    3) ~~~~~~~~~~~~~   final cull:  pCell->state.neighbors:" << std::endl;
//         for (const auto& nbr : pCell->state.neighbors)
//             std::cout << "             ID= " << pCell->ID<<" ("<<pCell->position[0]<<","<<pCell->position[1]<< ")" << " has nbr->ID = "<< nbr->ID <<" ("<<nbr->position[0]<<","<<nbr->position[1]<< ")" << std::endl;
//     }

//     pCell->state.neighbors.erase
//     (
//     std::remove_if(
//         pCell->state.neighbors.begin(),
//         pCell->state.neighbors.end(),
//         [&](const Cell* nbr) {          // <-- replace CellType* with your actual type
//             // double rest_length = (pCell->radius + nbr->radius) * rest_length_factor;
//             double rest_length = (pCell->phenotype.geometry.radius + nbr->phenotype.geometry.radius) * rest_length_factor;
//             double dx = nbr->position[0] - pCell->position[0];
//             double dy = nbr->position[1] - pCell->position[1];
//             double distance = std::sqrt(dx*dx + dy*dy);
//             return distance > rest_length;
//         }),
//     pCell->state.neighbors.end()
//     );


//     if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
//     {
//         std::cout << "    after 3) ~~~~~~~~~~~~~   final cull:  pCell->state.neighbors:" << std::endl;
//         for (const auto& nbr : pCell->state.neighbors)
//             std::cout << "             ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
//     }

//     // if (pCell->ID == 1 && PhysiCell_globals.current_time >= tmin && PhysiCell_globals.current_time <= tmax )
//     // {
//     //     std::cout << " ----- ID 1 has true nbrs, overlapping :" << std::endl;
//     //     for (const auto& nbr : pCell->state.neighbors)
//     //         std::cout <<  nbr->ID << std::endl;
//     // }
// }

//-----------------------------------
// void standard_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
// {
// 	if( pCell->functions.add_cell_basement_membrane_interactions )
// 	{
// 		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
// 	}
	
// 	pCell->state.simple_pressure = 0.0; 
// 	pCell->state.neighbors.clear(); // new 1.8.0
	
// 	//------------  First check the neighbors in my current voxel
// 	std::vector<Cell*>::iterator neighbor;
// 	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
// 	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
// 	{
// 		pCell->add_potentials(*neighbor);
// 	}
// 	std::vector<int>::iterator neighbor_voxel_index;
// 	std::vector<int>::iterator neighbor_voxel_index_end = 
// 		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();


// 	//------------  Then check the neighbors in surrounding voxels

// 	for( neighbor_voxel_index = 
// 		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
// 		neighbor_voxel_index != neighbor_voxel_index_end; 
// 		++neighbor_voxel_index )
// 	{
// 		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
// 			continue;

// 		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();

// 		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
// 		{
// 			pCell->add_potentials(*neighbor);
// 		}
// 	}

// 	pCell->update_motility_vector(dt); 
// 	pCell->velocity += phenotype.motility.motility_vector; 
	
// 	return; 
// }
//----------------------------------

// do every mechanics dt
void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    pCell->custom_data["cell_ID"] = pCell->ID;
    // static double effective_repulsion = 10.0;
    static double effective_repulsion = phenotype.mechanics.cell_cell_repulsion_strength;  // e.g., 10

    // generate_neighbor_list(pCell);   // <---------- NEW (rwh) - comment out and try using "nearby_cells" instead
    // double tmin=4;
    // double tmax= 6;
    // pCell->nearby_interacting_cells();
    // for (const auto& nbr : pCell->state.neighbors)

    if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    {
        std::cout << "------ nearby_cells\n";
        for (const auto& nbr : pCell->nearby_cells ())   // nearby_cells or nearby_interacting_cells
            std::cout << __FUNCTION__ << ": t=" <<PhysiCell_globals.current_time <<",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;

        std::cout << "------ nearby_interacting_cells\n";
        for (const auto& nbr : pCell->nearby_interacting_cells())   
            std::cout << __FUNCTION__ << ": t=" <<PhysiCell_globals.current_time <<",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
        // std::cout << " ~~nearby_interacting_cells are " << pCell->nearby_interacting_cells() << std::endl;
    }

	pCell->state.neighbors.clear();
    for (const auto& nbr : pCell->nearby_interacting_cells())  // nearby_cells or nearby_interacting_cells
    {
        if (pCell->ID == nbr->ID)
            continue;
        double rest_length = (pCell->phenotype.geometry.radius + nbr->phenotype.geometry.radius) * rest_length_factor;
        double dx = nbr->position[0] - pCell->position[0];
        double dy = nbr->position[1] - pCell->position[1];
        double distance = std::sqrt(dx*dx + dy*dy);
        if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
        {
            std::cout << "   t=" <<PhysiCell_globals.current_time <<" :  " << pCell->ID << " ->"<< nbr->ID << ": dist=" << distance <<" , rest_length="<<rest_length << std::endl;
        }
        if (distance <= rest_length)
        {
            pCell->state.neighbors.push_back(nbr);
        }
    }
    if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    {
        std::cout << "          -- state.neighbors  after constraining" << std::endl;
        for (const auto& nbr : pCell->state.neighbors)
            std::cout << ",         ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
    }

    // if (pCell->ID == 1 && PhysiCell_globals.current_time > tmin && PhysiCell_globals.current_time < tmax )
    // {
    //     std::cout << "    after 3) ~~~~~~~~~~~~~   final cull:  pCell->state.neighbors:" << std::endl;
    //     for (const auto& nbr : pCell->state.neighbors)
    //         std::cout << "             ID= " << pCell->ID << " has nbr->ID = " << nbr->ID << std::endl;
    // }


    // Chaste RepulsionForce / GeneralisedLinearSpringForce cell-cell repulsion.
    // Parameters matching Relax11.cpp:
    //   mu    = 5.0  (SetMeinekeSpringStiffness)
    //   alpha = 5.0  (hard-coded in GeneralisedLinearSpringForce.cpp, unused in pure repulsion)
    // Note: Chaste uses damping = 1.0 by default, so velocity = force / 1.0 = force.
    //       Set PhysiCell's damping coefficient to 1.0 in the XML config to match.
    // static double mu    = 5.0;
    // static double alpha = 5.0;  // exponential decay constant (attraction branch only)

    // static double mStiffness = 10;  // 5.0;
    // static double m_stiffness = parameters.doubles("m_stiffness");   // 5, 10, ?

    double r_a = phenotype.geometry.radius;

    // Reset velocity; motility not included here (matches Chaste's RelaxationForce-only setup)
    pCell->velocity = {0.0, 0.0, 0.0};

    pCell->custom_data["num_nbrs"] = 0;
    // int nbr_ID = -1;
    // double t_saved;
    // for( auto pNeighbor : *PhysiCell::all_cells )   // rwh: optimize
    for (const auto& pNeighbor : pCell->state.neighbors)
    {
        // if( pNeighbor == pCell ) continue;

        double r_b = pNeighbor->phenotype.geometry.radius;
        double rest_length = r_a + r_b;

        // if(std::abs(pCell->position[0]-pNeighbor->position[0]) > rest_length) continue;
        // else if(std::abs(pCell->position[1]-pNeighbor->position[1]) > rest_length) continue;

        // if (pCell->ID <= 1)
        // {
        //     std::cout << "---   nbr ID= " << pNeighbor->ID << ", rest_length= " << rest_length << std::endl;
        // }

        // Vector from A (pCell) toward B (pNeighbor)
        double dx = pNeighbor->position[0] - pCell->position[0];
        // if (PhysiCell_globals.current_time < 0.05 && (pCell->ID == 0) && (pNeighbor->ID==1))
        // {
        //     std::cout << "------ " << PhysiCell_globals.current_time << ", rest_length= " << rest_length << ", dx= "<<dx<< std::endl;   
        //     // ------ 0, rest_length= 10, dx= 5
        //     // ------ 0.01, rest_length= 10, dx= 4.99625
        //     // ------ 0.02, rest_length= 10, dx= 4.99374
        // }
        double dy = pNeighbor->position[1] - pCell->position[1];
        // double dy = 0.;
        double dz = 0.; // pNeighbor->position[2] - pCell->position[2];
        //jdouble distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        // double distance = std::sqrt( dx*dx + dy*dy);
        double distance = std::abs(dx);

        if( distance < 1e-12 ) continue;  // avoid division by zero

        double overlap = distance - rest_length;  // negative when cells overlap

        if( overlap >= 0.0 )
        {
            continue;

            // Attraction branch (linear * exponential decay); only reached if
            // distance >= rest_length -- RepulsionForce skips this, but kept
            // for completeness to match GeneralisedLinearSpringForce exactly.

            // double magnitude = mu * overlap * std::exp( -alpha * overlap / rest_length );
            // pCell->velocity[0] += magnitude * dx / distance;
            // pCell->velocity[1] += magnitude * dy / distance;
            // pCell->velocity[2] += magnitude * dz / distance;
        }
        else
        {
            pCell->custom_data["num_nbrs"] += 1;   // now performed in custom_cell_rule[_slow], after update_pos

            // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
            // {
            //     // std::cout << "---   cell ID= " << pCell->ID << " has just 1 nbr: ID = " << pNeighbor->ID << ", t="<<PhysiCell_globals.current_time << std::endl;
            //     nbr_ID = pNeighbor->ID;
            //     t_saved = PhysiCell_globals.current_time;
            //     // true_nbrs.push_back();

            //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
            //     // if (pNeighbor->ID != 0)
            //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
            // }
            // Repulsion branch: overlap/rest_length in (-1, 0)
            // log argument in (0, 1) -> magnitude < 0 -> force points away from B

            // logarithmic (matches Chaste's "stable" form)
            // double magnitude = mu * std::log( 1.0 + overlap / rest_length );

            // quadratic magnitude, preserving sign for direction
            // double magnitude = m_stiffness * overlap * std::abs(overlap) / rest_length;

            //------- PhysiCell
            // double temp_r = -distance; // -d

            // double temp_r = -distance; // -d
            // temp_r /= rest_length; // -d/R
            // temp_r += 1.0; // 1-d/R
            // temp_r *= temp_r; // (1-d/R)^2 

            // double effective_repulsion = sqrt( phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 
            // static double effective_repulsion = phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength ); 

            // temp_r *= effective_repulsion; 


            double temp = 1.0 - (distance / rest_length);   // normalized overlap fraction
            // double magnitude = stiffness * temp * temp;
            double magnitude = effective_repulsion * temp * temp;


            // --- cf. PhysiCell /core: void Cell::add_potentials(Cell* other_agent)
            //	velocity[i] += displacement[i] * temp_r; 
    	    // axpy( &velocity , temp_r , displacement );    # core/PhysiCell_cell.cpp

            // quadratic repulsion
            // pCell->velocity[0] += magnitude * dx / distance;
            // pCell->velocity[1] += magnitude * dy / distance;
            // pCell->velocity[2] += magnitude * dz / distance;
            // pCell->velocity[2] = 0.0;

            // PhysiCell
            // temp_r /= distance;
            // for( int i = 0 ; i < 3 ; i++ ) 
            // {
            //	velocity[i] += displacement[i] * temp_r; 
            // }
            // axpy( &velocity , temp_r , displacement ); 
            // pCell->velocity[0] += temp_r * dx / distance;
            // pCell->velocity[1] += temp_r * dy / distance;
            // pCell->velocity[0] += temp_r * dx;
            // pCell->velocity[0] += magnitude * dx / distance;
            pCell->velocity[0] -= magnitude * dx / distance;     // rwh: yipeee: negate for repulsion

            // pCell->velocity[1] += temp_r * dy;
            // pCell->velocity[1] = 0.0;
        }
    }
    // if ((pCell->custom_data["num_nbrs"] == 1) && pCell->ID == 1)
    // {
    //     std::cout << "--- FINAL:   cell ID= " << pCell->ID << " has 1 nbr: ID = " << nbr_ID << ", t_saved="<< t_saved << std::endl;
    //     // std::cout << "---   cell ID= " << pCell->ID << " had nbrs = " << pCell->state.neighbors <<", latest has ID = " << pNeighbor->ID << std::endl;
    //     // if (pNeighbor->ID != 0)
    //         // std::cout << "---   cell ID= " << pCell->ID << " has nbr ID = " << pNeighbor->ID << std::endl;
    // }
}

// called every dt_mech; pC->functions.custom_cell_rule( pC,pC->phenotype,time_since_last_mechanics );
void custom_cell_rule( Cell* pCell, Phenotype& phenotype , double dt )
{ 
    static bool reached_90 = false;
    std::stringstream ss;

    pCell->custom_data["cell_ID"] = pCell->ID;
    pCell->custom_data["num_nbrs"] = pCell->state.neighbors.size();
    pCell->custom_data["vel_mag"] = std::sqrt( pCell->previous_velocity[0]*pCell->previous_velocity[0] + pCell->previous_velocity[1]*pCell->previous_velocity[1] );


    // if (pCell->ID == 10 && pCell->position[0] >= 90.0)  // 90%, 90pct
    if (pCell->ID == 10)
    {
        // std::cout << "------- " << __FUNCTION__ << ":  ID= " << pCell->ID <<":  x= " << pCell->position[0] << std::endl;
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // full width: -5*d to 5*d, e.g. d=10 -> 50-(-50)=100; 90% is 90; divide by 2 = 45 (x-position)
        // if (pCell->position[0] >= 45.0)  // for symmetric test (hard-coded for diam=10)
        if (pCell->position[0] >= xpos_90pct)  // for symmetric test (hard-coded for diam=10)
        {
            if (!reached_90)
            {
                std::cout <<"---- "<< __FUNCTION__ << ": cell radius = " << phenotype.geometry.radius << std::endl;
                time_90pct = PhysiCell_globals.current_time;
                std::cout <<"---- "<< __FUNCTION__ << ": Width reached 90% , t= " << time_90pct << std::endl;

                std::ofstream outFile;
                // sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );
                ss << PhysiCell_settings.folder.c_str()  << "/time_90pct.txt" ;
                std::string out_file = ss.str();
                std::cout <<"---- "<< __FUNCTION__ << ": time_pct out_file=" << out_file << std::endl;
                outFile.open(out_file);   // read by  ../analysis/plot_11cells_sweep.py  for plotting results
                outFile << time_90pct;
                outFile.close();

                reached_90 = true;
                // std::exit();
            }
            else
            {
                // if (PhysiCell_globals.current_time / time_90pct > 10.0)
                if (cells11_flag && PhysiCell_globals.current_time / time_90pct > 10.0)
                {
                    std::cout <<"---- "<< __FUNCTION__ << "  calibrated T > 10.  Exit simulation!" << std::endl;
                    std::exit(-1);
                }

            }
        }
    }
} 