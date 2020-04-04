/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "memory_mgt.h"
#include "global_data.h"
#include "cor_global_data.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


Field_mem_info::Field_mem_info(const char *field_name, int decomp_id, int comp_or_grid_id, 
                               int buf_mark, const char *unit, const char *data_type, const char *annotation, bool check_field_name)
{
    int mem_size = 0;
    Remap_grid_class *remap_grid_grid = NULL, *remap_grid_decomp = NULL;
    Remap_data_field *remap_data_field;
	Decomp_info *decomp_info = NULL;

	
	total_dim_size_before_H2D = 1;
	total_dim_size_after_H2D = 1;

	chunk_field_instance_size = 0;

    if (decomp_id == -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, data_type != NULL, "Software error in Field_mem_info::Field_mem_info");
        if (comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,false)) {
            comp_id = comp_or_grid_id;
            grid_id = -1;
	        mem_size = get_data_type_size(data_type);
        }
        else {
            comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);
            grid_id = comp_or_grid_id;
            //EXECUTION_REPORT(REPORT_ERROR, comp_id, original_grid_mgr->get_original_CoR_grid(grid_id)->get_num_dimensions() == 1 && original_grid_mgr->get_original_CoR_grid(grid_id)->has_grid_coord_label(COORD_LABEL_LEV), 
            //                 "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": when the given parallel decomposition ID (the parameter \"decomp_id\") is -1, the corresponding grid \"%s\" must be an one-dimension vertical grid. Please check the model code with the annotation \"%s\"", field_name, original_grid_mgr->get_original_grid(comp_or_grid_id)->get_grid_name(), annotation);
            mem_size = original_grid_mgr->get_grid_size(grid_id, "in Field_mem_info::Field_mem_info") * get_data_type_size(data_type);
			remap_grid_grid = original_grid_mgr->get_original_CoR_grid(comp_or_grid_id);
        }
        host_comp_id = comp_id;
    }
    else {
		decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error2 in new Field_mem_info");
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Software error3 in new Field_mem_info");
        grid_id = comp_or_grid_id;
        comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);
        host_comp_id = decomp_info->get_host_comp_id();
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == decomps_info_mgr->get_comp_id_of_decomp(decomp_id), 
                         "Software error4 in new Field_mem_info");
        remap_grid_decomp = decomps_info_mgr->get_CoR_grid_of_decomp(decomp_id);
        remap_grid_grid = original_grid_mgr->get_original_CoR_grid(comp_or_grid_id);
		if (data_type != NULL && decomp_info->get_num_chunks() == 0)
	        mem_size = decomp_info->get_num_local_cells() * get_data_type_size(data_type) * remap_grid_grid->get_grid_size()/remap_grid_decomp->get_grid_size();
        EXECUTION_REPORT(REPORT_ERROR, host_comp_id, remap_grid_decomp->is_subset_of_grid(remap_grid_grid), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the grid (\"%s\") corresponding to the parallel decomposition \"decomp_id\" is not a subgrid of the grid (\"%s\") corresponding \"comp_or_grid_id\. Please check the model code with the annotation \"%s\"", field_name, decomp_info->get_grid_name(), original_grid_mgr->get_original_grid(comp_or_grid_id)->get_grid_name(), annotation);
    }
    
    host_comp_time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);

    const field_attr *field_attributes = fields_info->search_field_info(field_name);

    if (check_field_name) {
        EXECUTION_REPORT(REPORT_ERROR, host_comp_id, field_attributes != NULL, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the field name \"%s\" is unknown (has not been registered through the configuration XML file public_field_attribute.xml). Please check the model code with the annotation \"%s\"", field_name, field_name, annotation);
        bool dimensions_match_grid;
        if (words_are_the_same(field_attributes->field_dim, FIELD_0_DIM))
            dimensions_match_grid = decomp_id == -1 && remap_grid_grid == NULL;
        if (words_are_the_same(field_attributes->field_dim, FIELD_V1_DIM))
            dimensions_match_grid = decomp_id == -1 && remap_grid_grid != NULL && remap_grid_grid->get_num_dimensions() == 1 && remap_grid_grid->has_grid_coord_label(COORD_LABEL_LEV);
        else if (words_are_the_same(field_attributes->field_dim, FIELD_2_DIM))
            dimensions_match_grid = decomp_id != -1 && remap_grid_grid != NULL && remap_grid_grid->get_is_sphere_grid();
        else if (words_are_the_same(field_attributes->field_dim, FIELD_3_DIM))
            dimensions_match_grid = decomp_id != -1 && remap_grid_grid != NULL && remap_grid_grid->get_num_dimensions() == 3;
        else if (words_are_the_same(field_attributes->field_dim, FIELD_4_DIM))
            dimensions_match_grid = decomp_id != -1 && remap_grid_grid != NULL && remap_grid_grid->get_num_dimensions() == 4;
        if (!dimensions_match_grid) {
            if (grid_id != -1)
                EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\" at the model code with the annotation \"%s\": the dimension information (\"%s\") of the field that is specified in a configuration file does not match the dimensions of the corresponding grid \"%s\"", field_name, annotation, field_attributes->field_dim, original_grid_mgr->get_name_of_grid(grid_id));
            else EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the dimension information (\"%s\") of the field that is specified in a configuration file does not match the dimensions of the corresponding empty grid", field_name, annotation, field_attributes->field_dim);
        }
    }    

    if (strlen(unit) > 0)
        strcpy(field_unit, unit);
    else if (field_attributes != NULL) 
        strcpy(field_unit, fields_info->search_field_info(field_name)->field_unit);
    // check the field unit

    strcpy(this->field_name, field_name);
    this->decomp_id = decomp_id;
    this->comp_or_grid_id = comp_or_grid_id;
    this->buf_mark = buf_mark;
    this->usage_tag = -1;
    is_registered_model_buf = false;
    is_field_active = false;
    define_order_count = -1;
    last_define_time = 0x7fffffffffffffff;

    remap_data_field = new Remap_data_field;
    strcpy(remap_data_field->field_name_in_application, field_name);
    strcpy(remap_data_field->field_name_in_IO_file, field_name);
	num_chunks = 0;
	chunks_buf = NULL;
	chunks_data_buf_size = NULL;
	if (decomp_info != NULL && decomp_info->get_num_chunks() > 0) {
		chunks_buf = new void* [decomp_info->get_num_chunks()];
		chunks_data_buf_size = new int [decomp_info->get_num_chunks()];
		for (int i = 0; i < decomp_info->get_num_chunks(); i ++)
			chunks_buf[i] = NULL;
		if (data_type != NULL) {
			num_chunks = decomp_info->get_num_chunks();
			for (int i = 0; i < num_chunks; i ++) {
				chunks_data_buf_size[i] = decomp_info->get_chunk_size(i)*remap_grid_grid->get_grid_size()/remap_grid_decomp->get_grid_size();
				if (decomp_info->get_chunk_size(i) > 0) 
					chunks_buf[i] = new char [chunks_data_buf_size[i]*get_data_type_size(data_type)];
				chunk_field_instance_size += chunks_data_buf_size[i]; 
			}
		}
	}
	if (data_type != NULL) {
	    strcpy(remap_data_field->data_type_in_application, data_type);
	    strcpy(remap_data_field->data_type_in_IO_file, data_type);
	    remap_data_field->required_data_size = mem_size / get_data_type_size(data_type);
	    if (remap_data_field->required_data_size > 0) {
	        remap_data_field->data_buf = (char*) (new long [(mem_size+sizeof(long)-1)/sizeof(long)]);
	        memset(remap_data_field->data_buf, 0, mem_size);
	    }
		if (chunk_field_instance_size > 0)
			remap_data_field->required_data_size = chunk_field_instance_size;
	    remap_data_field->read_data_size = remap_data_field->required_data_size;
	}
    if (check_field_name)
        remap_data_field->set_field_long_name(fields_info->get_field_long_name(field_name));
    remap_data_field->set_field_unit(unit);   // to complete: when strlen(unit) is 0, use default unit of the field

    if (grid_id == -1)
        grided_field_data = new Remap_grid_data_class(NULL, remap_data_field);
	else if (decomp_id == -1)
		grided_field_data = new Remap_grid_data_class(remap_grid_grid, remap_data_field);
    else {
        Remap_grid_class *decomp_grid = decomp_grids_mgr->search_decomp_grid_info(decomp_id, remap_grid_grid, false)->get_decomp_grid();
        grided_field_data = new Remap_grid_data_class(decomp_grid, remap_data_field);
		if (data_type != NULL)
	        remap_data_field->set_fill_value(NULL);
    }
    
    last_checksum = -1;

	if (grid_id != -1)
		remap_grid_grid->get_total_dim_size_before_and_after_H2D(total_dim_size_before_H2D, total_dim_size_after_H2D);
	if (words_are_the_same(field_name, V3D_GRID_3D_LEVEL_FIELD_NAME))
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grided_field_data->get_coord_value_grid()->get_num_dimensions() == 3 && grided_field_data->get_coord_value_grid()->does_use_V3D_level_coord());
}


Field_mem_info::~Field_mem_info()
{
    if (is_registered_model_buf && num_chunks == 0)
        grided_field_data->get_grid_data_field()->data_buf = NULL;
    delete grided_field_data;
	if (chunks_buf != NULL) {
		if (!is_registered_model_buf)
			for (int i = 0; i < num_chunks; i ++)
				if (chunks_buf[i] != NULL)
					delete [] chunks_buf[i];
		delete [] chunks_buf;
		delete [] chunks_data_buf_size;
	}
}


void Field_mem_info::reset_mem_buf(void * buf, bool is_external_field, int usage_tag)
{
	if (get_size_of_field() > 0) {
	    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, buf != NULL, "The data buffer corresponding to the field instance of \"%s\" is not allocated. Please verify the model code corresponding to the annotation \"%s\"", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
		
	}
    EXECUTION_REPORT(REPORT_ERROR, -1, !is_registered_model_buf, "Software error to release a registered buffer");

    if (grided_field_data->get_grid_data_field()->data_buf != NULL)
        delete [] grided_field_data->get_grid_data_field()->data_buf;

    grided_field_data->get_grid_data_field()->data_buf = buf;

    if (is_external_field) {
        this->usage_tag = usage_tag;
		is_registered_model_buf = true;
    }
}


void Field_mem_info::reset_mem_buf(Field_mem_info *original_field_mem, int usage_tag)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !is_registered_model_buf, "Software error to release a registered buffer");
	if (original_field_mem->num_chunks > 0)
	    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_field_mem->grided_field_data->get_grid_data_field()->data_buf == NULL, "Software error to release a registered buffer");

	if (grided_field_data->get_grid_data_field()->data_buf != NULL) {
		delete [] grided_field_data->get_grid_data_field()->data_buf;
		grided_field_data->get_grid_data_field()->data_buf = NULL;
	}
	this->num_chunks = original_field_mem->num_chunks;
	this->chunk_field_instance_size = original_field_mem->chunk_field_instance_size;
	if (this->num_chunks > 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->chunks_buf != NULL, "Software error in Field_mem_info::reset_mem_buf");
		for (int i = 0; i < num_chunks; i ++) {
			if (this->chunks_buf[i] != NULL)
				delete [] this->chunks_buf[i];
			this->chunks_buf[i] = original_field_mem->chunks_buf[i];
			this->chunks_data_buf_size[i] = original_field_mem->chunks_data_buf_size[i];
		}
	}
	else this->grided_field_data->get_grid_data_field()->data_buf = original_field_mem->grided_field_data->get_grid_data_field()->data_buf;
    this->usage_tag = usage_tag;
	this->is_registered_model_buf = original_field_mem->is_registered_model_buf;
}


void Field_mem_info::change_datatype_to_double()
{
    grided_field_data->change_datatype_in_application(DATA_TYPE_DOUBLE);
}


void Field_mem_info::define_field_values(bool is_restarting)
{
    if (!is_restarting)
        is_field_active = true;
    last_define_time = host_comp_time_mgr->get_current_full_time();
}


void Field_mem_info::use_field_values(const char *annotation)
{    
    if (is_registered_model_buf) 
        return;
    
    if (last_define_time == host_comp_time_mgr->get_current_full_time())
        return;

    if (last_define_time == 0x7fffffffffffffff) {
        if (is_registered_model_buf)
            EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "field instance (field_name=\"%s\", decomp_name=\"%s\", grid_name=\"%s\", bufmark=%x) is used before defining it. Please modify the model code with the annotation \"%s\"", field_name, get_decomp_name(), get_grid_name(), buf_mark, annotation);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Field_mem_info::use_field_values: field instance (field_name=\"%s\", decomp_name=\"%s\", grid_name=\"%s\", bufmark=%x) is used before defining it. Please modify the model code with the annotation \"%s\"", field_name, get_decomp_name(), get_grid_name(), buf_mark, annotation);        
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, last_define_time <= host_comp_time_mgr->get_current_full_time(), "C-Coupler error in Field_mem_info::use_field_values: wrong time order");
}


bool Field_mem_info::match_field_instance(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark)
{
    return words_are_the_same(this->field_name, field_name) && this->decomp_id == decomp_id && this->comp_or_grid_id == comp_or_grid_id && this->buf_mark == buf_mark;
}


bool Field_mem_info::match_field_mem(void *data_buffer)
{
    return this->get_data_buf() == data_buffer;
}


void Field_mem_info::reset_field_name(const char *new_name)
{
    strcpy(field_name, new_name);
}


void Field_mem_info::calculate_field_conservative_sum(Field_mem_info *area_field)
{
    double partial_sum, total_sum;
    long size;

    if (report_internal_log_enabled) {
        EXECUTION_REPORT(REPORT_ERROR,-1, words_are_the_same(get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "C-Coupler error in calculate_field_sum");
        size = get_field_data()->get_grid_data_field()->required_data_size;
        partial_sum = 0;
        for (long j = 0; j < size; j ++)
            partial_sum += (((double*) get_data_buf())[j])*(((double*) area_field->get_data_buf())[j]);
        MPI_Allreduce(&partial_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(area_field->get_comp_id(),"in Field_mem_info::calculate_field_conservative_sum"));
    }
}


void Field_mem_info::get_total_dim_size_before_and_after_H2D(int &total_dim_size_before_H2D, int &total_dim_size_after_H2D)
{
	total_dim_size_before_H2D = this->total_dim_size_before_H2D;
	total_dim_size_after_H2D = this->total_dim_size_after_H2D;	
}


void Field_mem_info::check_field_sum(bool do_check_sum, bool bypass_decomp, const char *hint)
{
    int total_dim_size_before_H2D = 1, total_dim_size_after_H2D = 1, decomp_size = 1;
	const int *decomp_local_cell_global_indx = NULL;
    long size;
	unsigned long partial_sum = 0, total_sum;
	char *current_data_buf;


	if (!do_check_sum)
		return;

	if (chunks_buf != NULL && !is_registered_model_buf)
		return;
	
    if (report_error_enabled && is_registered_model_buf) {
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
        EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "Try to check the model data buffers of the field \"%s\" registered corresponding to the code annotation \"%s\". If it fails to pass the check (the model run is stopped), please make sure corresponding model data buffers are global variables and have not been released", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
		if (get_size_of_field() > 0) {
	        char *temp_array = new char [get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application)];
			if (get_data_buf() != NULL) {
		        memcpy(temp_array, get_data_buf(), get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application));
		        memcpy(get_data_buf(), temp_array, get_size_of_field()*get_data_type_size(get_field_data()->get_grid_data_field()->data_type_in_application));
			}
			if (num_chunks > 0) {
				for (int i = 0; i < num_chunks; i ++) {
					Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
					if (chunks_buf[i] != NULL) {
						memcpy(temp_array, chunks_buf[i], chunk_field_instance_size/decomp_info->get_num_local_cells()*decomp_info->get_chunk_size(i)*get_data_type_size(get_data_type()));
						memcpy(chunks_buf[i], temp_array, chunk_field_instance_size/decomp_info->get_num_local_cells()*decomp_info->get_chunk_size(i)*get_data_type_size(get_data_type()));
					}
				}
			}
	        delete [] temp_array;
		}
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));        
        EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "Pass the check of the model data buffers of the field \"%s\" registered corresponding to the code annotation \"%s\". If it fails to pass the check (the model run is stopped), please make sure corresponding model data buffers are global variables and have not been released", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
    }

	if (grid_id != -1) {
		get_total_dim_size_before_and_after_H2D(total_dim_size_before_H2D, total_dim_size_after_H2D);
		if (decomp_id != -1) {
			decomp_size = decomps_info_mgr->get_decomp_info(decomp_id)->get_num_local_cells();
			decomp_local_cell_global_indx = decomps_info_mgr->get_decomp_info(decomp_id)->get_local_cell_global_indx();
		}
	}
	total_dim_size_before_H2D *= get_data_type_size(grided_field_data->get_grid_data_field()->data_type_in_application);

	if (num_chunks == 0) {
		for (int k = 0; k < total_dim_size_after_H2D; k ++) {	
			for (int j = 0; j < decomp_size; j ++) {
				if (!bypass_decomp && decomp_local_cell_global_indx != NULL && decomp_local_cell_global_indx[j] == CCPL_NULL_INT)
					continue;
				current_data_buf = (char*)get_data_buf() + k*decomp_size*total_dim_size_before_H2D+total_dim_size_before_H2D*j;
				for (int i = 0; i < total_dim_size_before_H2D; i ++)
					partial_sum += (((unsigned long)(current_data_buf[i])) << ((i%8)*8));
			}
		}
	}
	else {
		Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
		for (int m = 0; m < num_chunks; m ++) {
			for (int k = 0; k < total_dim_size_after_H2D; k ++) {	
				for (int j = 0; j < decomp_info->get_chunk_size(m); j ++) {
					if (!bypass_decomp && decomp_local_cell_global_indx != NULL && decomp_local_cell_global_indx[decomp_info->get_chunk_start(m)+j] == CCPL_NULL_INT)
						continue;
					current_data_buf = (char*)(chunks_buf[m]) + k*decomp_info->get_chunk_size(m)*total_dim_size_before_H2D+total_dim_size_before_H2D*j;
					for (int i = 0; i < total_dim_size_before_H2D; i ++)
						partial_sum += (((unsigned long)(current_data_buf[i])) << ((i%8)*8));
				}
			}
		}
	}

    if (decomp_id != -1)
        MPI_Allreduce(&partial_sum, &total_sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
	else {
		total_sum = partial_sum;
		MPI_Bcast(&total_sum, 1, MPI_UNSIGNED_LONG, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "Field_mem_info::check_field_sum"));
		EXECUTION_REPORT(REPORT_WARNING, host_comp_id, partial_sum == total_sum, "As an instance of the field \"%s\" is not on a horizontal grid, all its values should be the same but currently are not the same across all processes of the corresponding component model. Please check the model code related to the annotation \"%s\"", field_name, annotation_mgr->get_annotation(field_instance_id, "allocate field instance"));
	}
	if (comp_comm_group_mgt_mgr->search_global_node(host_comp_id) != NULL)
		if (comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_current_proc_local_id() == 0)
			EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "Check sum of field \"%s\" %s is %llx", get_field_name(), hint, total_sum);
}


bool Field_mem_info::field_has_been_defined()
{
    return last_define_time != 0x7fffffffffffffff;
}


long Field_mem_info::get_size_of_field()
{
	if (num_chunks > 0) 
		return chunk_field_instance_size;
	if (decomp_id == -1)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grided_field_data->get_grid_data_field()->required_data_size > 0, "Software error in Field_mem_info::get_size_of_field");
    return grided_field_data->get_grid_data_field()->required_data_size;
}


const char *Field_mem_info::get_grid_name() 
{
    if (grid_id == -1)
        return NULL;

    return original_grid_mgr->search_grid_info(grid_id)->get_grid_name();
}


const char *Field_mem_info::get_decomp_name()
{
    if (decomp_id == -1)
        return NULL;

    return decomps_info_mgr->get_decomp_info(decomp_id)->get_decomp_name();
}


const char *Field_mem_info::get_data_type()
{
    return get_field_data()->get_grid_data_field()->data_type_in_application;
}


void Field_mem_info::set_field_instance_id(int field_instance_id, const char *annotation)
{
    this->field_instance_id = field_instance_id;
    annotation_mgr->add_annotation(field_instance_id, "allocate field instance", annotation);
}


long Field_mem_info::calculate_overall_checksum()
{
	if (num_chunks == 0)
		return calculate_checksum_of_array((const char*)get_data_buf(), get_size_of_field(), get_data_type_size(get_data_type()), NULL, NULL);

	long overall_checksum = 0;
	Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
	for (int i = 0; i < num_chunks; i ++)
		overall_checksum += calculate_checksum_of_array((char*)(chunks_buf[i]), chunk_field_instance_size/decomp_info->get_num_local_cells()*decomp_info->get_chunk_size(i), get_data_type_size(get_data_type()), NULL, NULL);

	return overall_checksum;
}


bool Field_mem_info::is_checksum_changed()
{
    if (last_checksum == -1)
        return true;

    long current_checksum = calculate_overall_checksum();
    
    return current_checksum != last_checksum;
}


void Field_mem_info::reset_checksum()
{
    last_checksum = calculate_overall_checksum();
}


void Field_mem_info::finish_chunk_registration(const char *annotation)
{
    synchronize_comp_processes_for_API(comp_id, API_ID_FIELD_MGT_FINISH_CHUNK_FIELD_INST_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"Field_mem_info::finalize_chunk_registration"), "Finish registering a chunk field instance", annotation);
	EXECUTION_REPORT(REPORT_ERROR, comp_id, decomp_id != -1, "ERROR happens when calling the API \"CCPL_finish_chunk_field_instance_registration\" for a field instance of \"%s\": it is not a chunk field instance but a scalar field instance. Please verify the model code with the annotation \"%s\"", field_name, annotation);
	Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
	EXECUTION_REPORT(REPORT_ERROR, comp_id, decomp_info->get_num_chunks() > 0, "ERROR happens when calling the API \"CCPL_finish_chunk_field_instance_registration\" for a field instance of \"%s\": it is not a chunk field instance because it does not correspond to a chunk parallel decomposition (\"%s\"). Please verify the model code with the annotation \"%s\"", field_name, decomp_info->get_decomp_name(), annotation);
	EXECUTION_REPORT(REPORT_ERROR, comp_id, decomp_info->get_num_chunks() == num_chunks, "ERROR happens when calling the API \"CCPL_finish_chunk_field_instance_registration\" for a field instance of \"%s\": the data buffers of some chunks of this field instance has not been specified. Please verify the model code with the annotation \"%s\"", field_name, annotation);
	is_registered_model_buf = true;
}


void Field_mem_info::register_a_chunk(void *buf, int field_size, const char *data_type, const char *annotation)
{
	Remap_grid_class *remap_grid_grid = NULL, *remap_grid_decomp = NULL;

	
	EXECUTION_REPORT(REPORT_ERROR, comp_id, decomp_id != -1, "ERROR happens when calling the API \"CCPL_register_chunk_of_field_instance\" for a field instance of \"%s\": it is not a chunk field instance but a scalar field instance. Please verify the model code with the annotation \"%s\"", field_name, annotation);
	Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
	EXECUTION_REPORT(REPORT_ERROR, comp_id, decomp_info->get_num_chunks() > 0, "ERROR happens when calling the API \"CCCPL_register_chunk_of_field_instance\" for a field instance of \"%s\": it is not a chunk field instance because it does not correspond to a chunk parallel decomposition (\"%s\"). Please verify the model code with the annotation \"%s\".", field_name, decomp_info->get_decomp_name(), annotation);
	EXECUTION_REPORT(REPORT_ERROR, comp_id, num_chunks < decomp_info->get_num_chunks(), "ERROR happens when calling the API \"CCCPL_register_chunk_of_field_instance\" for a field instance of \"%s\": the number of registered chunks is larger than the maximum number (%d). Please verify the model code with the annotation \"%s\".", field_name, decomp_info->get_num_chunks(), annotation);
	remap_grid_grid = original_grid_mgr->get_original_CoR_grid(grid_id);
	remap_grid_decomp = decomps_info_mgr->get_CoR_grid_of_decomp(decomp_id);
	if (decomp_info->get_chunk_size(num_chunks) > 0)
		EXECUTION_REPORT(REPORT_ERROR, comp_id, field_size == decomp_info->get_chunk_size(num_chunks)*remap_grid_grid->get_grid_size()/remap_grid_decomp->get_grid_size(), "ERROR happens when calling the API \"CCCPL_register_chunk_of_field_instance\" for registering No. %d chunk of a field instance of \"%s\": the size of data buffer (%d) is different from the required size (%d). Please verify the model code with the annotation \"%s\".", num_chunks+1, field_name, field_size, decomp_info->get_chunk_size(num_chunks)*remap_grid_grid->get_grid_size()/remap_grid_decomp->get_grid_size(), annotation);
	if (num_chunks == 0) {
	    strcpy(grided_field_data->get_grid_data_field()->data_type_in_application, data_type);
	    strcpy(grided_field_data->get_grid_data_field()->data_type_in_IO_file, data_type);
	}
	else EXECUTION_REPORT(REPORT_ERROR, comp_id, words_are_the_same(grided_field_data->get_grid_data_field()->data_type_in_application, data_type), "ERROR happens when calling the API \"CCCPL_register_chunk_of_field_instance\" for a field instance of \"%s\": the data types among chunks are not the same (\"%s\" vs \"%s\"). Please verify the model code with the annotation \"%s\".", field_name, grided_field_data->get_grid_data_field()->data_type_in_application, data_type, annotation);
	for (int i = 0; i < num_chunks; i ++)
		if (chunks_buf[i] != NULL || buf != NULL)
			EXECUTION_REPORT(REPORT_ERROR, comp_id, chunks_buf[i] != buf, "ERROR happens when calling the API \"CCCPL_register_chunk_of_field_instance\" for a field instance of \"%s\": the No. %d and No. %d chunks correspond to the same data buffer. Please verify the model code with the annotation \"%s\".", field_name, i+1, num_chunks+1, annotation);
	chunks_data_buf_size[num_chunks] = field_size;
	chunk_field_instance_size += field_size;
	if (decomp_info->get_chunk_size(num_chunks) == 0)
		chunks_buf[num_chunks++] = NULL;
	else chunks_buf[num_chunks++] = buf;
}


void Field_mem_info::confirm_overall_data_buf_for_chunks()
{
	if (num_chunks > 0 && grided_field_data->get_grid_data_field()->data_buf == NULL && chunk_field_instance_size > 0)
		grided_field_data->get_grid_data_field()->data_buf = new char [chunk_field_instance_size*get_data_type_size(get_data_type())];	
}


void Field_mem_info::transformation_between_chunks_array(bool chunks_to_array)
{
	if (num_chunks == 0)
		return;
	
	Decomp_info *decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
	long offset = 0;
	int data_type_size = get_data_type_size(get_data_type());

	confirm_overall_data_buf_for_chunks();
		
	for (int j = 0; j < total_dim_size_after_H2D; j ++) {
		for (int i = 0; i < num_chunks; i ++) {
			if (chunks_to_array)
				memcpy((char*)get_data_buf()+offset, (char*)(chunks_buf[i])+data_type_size*total_dim_size_before_H2D*decomp_info->get_chunk_size(i)*j, data_type_size*total_dim_size_before_H2D*decomp_info->get_chunk_size(i));
			else memcpy((char*)(chunks_buf[i])+data_type_size*total_dim_size_before_H2D*decomp_info->get_chunk_size(i)*j, (char*)get_data_buf()+offset, data_type_size*total_dim_size_before_H2D*decomp_info->get_chunk_size(i));
			offset += data_type_size*total_dim_size_before_H2D*decomp_info->get_chunk_size(i);
		}
	}
}


void Field_mem_info::change_to_registered_without_data_buffers()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !is_registered_model_buf, "Software error in Field_mem_info::change_to_registered_without_data_buffers");
	for (int i = 0; i < num_chunks; i ++)
		if (chunks_buf[i] != NULL) {
			delete [] chunks_buf[i];
			chunks_buf[i] = NULL;
		}
		
	if (grided_field_data->get_grid_data_field()->data_buf != NULL) {
		delete [] grided_field_data->get_grid_data_field()->data_buf;
		grided_field_data->get_grid_data_field()->data_buf = NULL;
	}
	
	is_registered_model_buf = true;
}


void Field_mem_info::set_data_buf_from_model(void *data_buf, int chunk_index)
{
	if (num_chunks > 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_registered_model_buf && chunk_index >= 0 && chunk_index <= num_chunks && chunks_buf[chunk_index] == NULL, "Software error in Field_mem_info::change_to_registered_without_data_buffers");
	}
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_registered_model_buf && chunk_index == -1 && grided_field_data->get_grid_data_field()->data_buf == NULL, "Software error in Field_mem_info::change_to_registered_without_data_buffers");
	if (chunk_index == -1)
		grided_field_data->get_grid_data_field()->data_buf = data_buf;
	else chunks_buf[chunk_index] = data_buf;
}



Field_mem_info *Memory_mgt::alloc_mem(Field_mem_info *original_field_instance, int special_buf_mark, int object_id, const char *unit_or_datatype, bool check_field_name)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, special_buf_mark == BUF_MARK_DATATYPE_TRANS || special_buf_mark == BUF_MARK_AVERAGED_INNER || special_buf_mark == BUF_MARK_AVERAGED_INTER || special_buf_mark == BUF_MARK_UNIT_TRANS || special_buf_mark == BUF_MARK_DATA_TRANSFER || 
                     special_buf_mark == BUF_MARK_IO_FIELD_MIRROR || special_buf_mark == BUF_MARK_REMAP_NORMAL || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST || special_buf_mark == BUF_MARK_REMAP_FRAC || special_buf_mark == BUF_MARK_AVERAGED_BACKUP || special_buf_mark == BUF_MARK_ENS_DATA_TRANSFER, "Software error in Field_mem_info *alloc_mem: wrong special_buf_mark");
    int new_buf_mark = (special_buf_mark ^ object_id);
    Field_mem_info *existing_field_instance = search_field_instance(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark);
    if (existing_field_instance != NULL) {
        if (special_buf_mark == BUF_MARK_UNIT_TRANS)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), original_field_instance->get_data_type()) && words_are_the_same(existing_field_instance->get_unit(), unit_or_datatype), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with different data type or wrong unit", new_buf_mark);
        else if (unit_or_datatype != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), unit_or_datatype), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with wrong data type", new_buf_mark);
        else EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(existing_field_instance->get_data_type(), original_field_instance->get_data_type()), "Software error in Field_mem_info *alloc_mem: special field instance exists %lx with wrong data type", new_buf_mark);
        return existing_field_instance;
    }
    if (special_buf_mark == BUF_MARK_AVERAGED_INNER || special_buf_mark == BUF_MARK_AVERAGED_INTER || special_buf_mark == BUF_MARK_AVERAGED_BACKUP)
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), original_field_instance->get_data_type(), "new field instance for averaging", check_field_name));    
    else if (special_buf_mark == BUF_MARK_REMAP_FRAC)
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), original_field_instance->get_data_type(), "new field instance for the remapping with fraction", check_field_name));
    else if (special_buf_mark == BUF_MARK_DATATYPE_TRANS || special_buf_mark == BUF_MARK_DATA_TRANSFER || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC || special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST) {
        get_data_type_size(unit_or_datatype);
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_IO_FIELD_MIRROR) {
        get_data_type_size(unit_or_datatype);
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation", check_field_name));        
    }
    else if (special_buf_mark == BUF_MARK_UNIT_TRANS) {
        // check unit


        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, unit_or_datatype, original_field_instance->get_data_type(), "new field instance for unit transformation", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_NORMAL) {
        // check unit


        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for remapping", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_SRC) {
        // check unit
    

    fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation in remapping", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_REMAP_DATATYPE_TRANS_DST) {
        // check unit


        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for data type transformation in remapping", check_field_name));
    }
    else if (special_buf_mark == BUF_MARK_ENS_DATA_TRANSFER) {
        fields_mem.push_back(new Field_mem_info(original_field_instance->get_field_name(), original_field_instance->get_decomp_id(), original_field_instance->get_comp_or_grid_id(), new_buf_mark, original_field_instance->get_unit(), unit_or_datatype, "new field instance for statistical processing in ensemble procedures", check_field_name));
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in Field_mem_info *alloc_mem");

    fields_mem[fields_mem.size()-1]->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|(fields_mem.size()-1), "in Memory_mgt::alloc_mem");
	

    return fields_mem[fields_mem.size()-1];
}


Field_mem_info *Memory_mgt::alloc_mem(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark, const char *data_type, const char *field_unit, const char *annotation, bool check_field_name)
{
    Field_mem_info *field_mem, *pair_field;
    int i, comp_id;
    bool find_field_in_cfg;


    EXECUTION_REPORT(REPORT_ERROR, -1, buf_mark < 0, "Software error in Memory_mgt::alloc_mem: wrong value of buffer mark");
    EXECUTION_REPORT(REPORT_ERROR, -1, data_type != NULL, "Software error in Memory_mgt::alloc_mem: data type is NULL");
    get_data_type_size(data_type);
    if (decomp_id != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(decomp_id), "Software error in Memory_mgt::alloc_mem: wrong decomposition id");
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error in Memory_mgt::alloc_mem: wrong grid id");
        comp_id = original_grid_mgr->search_grid_info(comp_or_grid_id)->get_comp_id();
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,false) || original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Software error in Memory_mgt::alloc_mem: wrong component id or grid_id");
		if (original_grid_mgr->is_grid_id_legal(comp_or_grid_id))
			comp_id = original_grid_mgr->search_grid_info(comp_or_grid_id)->get_comp_id();
		else comp_id = comp_or_grid_id;
    }

    
    /* If memory buffer has been allocated, return it */
    for (i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->match_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark)) {
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(data_type, fields_mem[i]->get_field_data()->get_grid_data_field()->data_type_in_application),
                             "Software error in Memory_mgt::alloc_mem: data types conflict");
            return fields_mem[i];
        }

    /* Compute the size of the memory buffer and then allocate and return it */
    field_mem = new Field_mem_info(field_name, decomp_id, comp_or_grid_id, buf_mark, field_unit, data_type, annotation, check_field_name);
    field_mem->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|fields_mem.size(), annotation);
    fields_mem.push_back(field_mem);

    return field_mem;
}


Memory_mgt::~Memory_mgt()
{
    for (int i = 0; i < fields_mem.size(); i ++)
        delete fields_mem[i];
}


Field_mem_info *Memory_mgt::search_field_via_data_buf(const void *data_buf, bool diag)
{
    for (int i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->get_data_buf() == data_buf)
            return fields_mem[i];

    if (diag)
        EXECUTION_REPORT(REPORT_ERROR,-1, false, "C-Coupler error in search_field_via_data_buf\n");
    return NULL;
}


void Memory_mgt::check_sum_of_all_registered_fields(int comp_id, bool bypass_decomp, const char *hint)
{
    for (int i = 0; i < fields_mem.size(); i ++)
		if (fields_mem[i]->get_comp_id() == comp_id && fields_mem[i]->get_is_registered_model_buf() && (fields_mem[i]->get_decomp_id() == -1 || decomps_info_mgr->get_decomp_info(fields_mem[i]->get_decomp_id())->is_registered_decomp()))
	        fields_mem[i]->check_field_sum(true, bypass_decomp, hint);
}


int Memory_mgt::get_field_size(void *data_buf, const char *annotation)
{
    Field_mem_info *field = search_field_via_data_buf(data_buf, false);

    EXECUTION_REPORT(REPORT_ERROR,-1, field != NULL, "Detect a memory buffer that is not managed by C-Coupler. Please verify the model code according to annotation \"%s\"", annotation);

    return field->get_size_of_field();
}


Field_mem_info *Memory_mgt::search_field_instance(const char *field_name, int decomp_id, int comp_or_grid_id, int buf_mark)
{
    for (int i = 0; i < fields_mem.size(); i ++)
        if (fields_mem[i]->match_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark))
            return fields_mem[i];

    return NULL;
}


int Memory_mgt::register_external_field_instance(const char *field_name, void *data_buffer, int field_size, int decomp_id, int comp_or_grid_id, 
                                                 int buf_mark, int usage_tag, const char *unit, const char *data_type, const char *annotation)
{
    int comp_id, API_id;
    Field_mem_info *existing_field_instance_instance, *new_field_instance;
	MPI_Comm comm;


    if (buf_mark == BUF_MARK_IO_FIELD_REG)
        API_id = API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER;
    else if (data_type != NULL)
		API_id = API_ID_FIELD_MGT_REG_FIELD_INST;
	else API_id = API_ID_FIELD_MGT_START_CHUNK_FIELD_INST_REG;
	check_API_parameter_comp_or_grid(comp_or_grid_id, API_id, MPI_COMM_NULL, "registering a field instance or a I/O field", "comp_or_grid_ID", annotation);
    if (comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_or_grid_id,true))
        comp_id = comp_or_grid_id;
    else comp_id = original_grid_mgr->get_comp_id_of_grid(comp_or_grid_id);
    check_API_parameter_string_length(comp_id, API_id, CCPL_NAME_STR_LEN, field_name, "field_name", annotation);
    check_and_verify_name_format_of_string_for_API(comp_id, field_name, API_id, "the field instance", annotation);
	comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id,"C-Coupler code in register_external_field_instance for getting component management node");
	check_API_parameter_decomp(comp_id, API_id, comm, "registering a field instance or a I/O field", decomp_id, "decomp_id", annotation);

    if (decomp_id != -1) {
        EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(comp_or_grid_id), "Error happens when calling the API \"CCPL_register_field_instance/CCPL_start_chunk_field_instance_registration\" to register a field instance of \"%s\": the parameter of grid ID is wrong. Please check the model code with the annotation \"%s\".", field_name, annotation);        
        EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == decomps_info_mgr->get_comp_id_of_decomp(decomp_id), 
                         "Error happens when calling the API \"CCPL_register_field_instance/CCPL_start_chunk_field_instance_registration\" to register a field instance of \"%s\": the parameters of grid ID and decomposition ID do not match each other: they belong to different component models. Please check the model code with the annotation \"%s\".",
                         field_name, annotation);
		if (data_type != NULL) 
			EXECUTION_REPORT(REPORT_ERROR, comp_id, decomps_info_mgr->get_decomp_info(decomp_id)->get_num_chunks() == 0, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the API \"CCPL_start_chunk_field_instance_registration\" should be used because the corresponding parallel decomposition \"%s\" is a chunk decomposition. Please check the model code with the annotation \"%s\".", field_name, decomps_info_mgr->get_decomp_info(decomp_id)->get_decomp_name(), annotation);
    }
    synchronize_comp_processes_for_API(comp_id, API_id, comm, "registering a field instance or a I/O field", annotation);
    comp_comm_group_mgt_mgr->confirm_coupling_configuration_active(comp_id, API_id, true, annotation);
    check_API_parameter_string(comp_id, API_id, comm, "registering a field instance or a I/O field", field_name, "field name", annotation);
    check_API_parameter_string(comp_id, API_id, comm, "registering a field instance or a I/O field", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id,true,"C-Coupler code in register_external_field_instance for getting component management node")->get_comp_name(), "the component name specified by the corresponding ID", annotation);
    check_API_parameter_int(comp_id, API_id, comm, NULL, buf_mark, "buf_mark", annotation);
    check_API_parameter_int(comp_id, API_id, comm, NULL, usage_tag, "usage_tag", annotation);
	if (data_type != NULL)
	    check_API_parameter_string(comp_id, API_id, comm, "registering a field instance or a I/O field", data_type, "the implicit data type (such as integer, float, and double) of the field instance", annotation);
    check_API_parameter_string(comp_id, API_id, comm, "registering a field instance or a I/O field", unit, "unit", annotation);
    if (decomp_id == -1)
        //EXECUTION_REPORT(REPORT_ERROR, -1, comp_id == comp_or_grid_id, "Error happens when calling the API \"CCPL_register_field_instance/CCPL_start_chunk_field_instance_registration\" to register a field instance of \"%s\" at the model code with the annotation \"%s\". We are sorry that C-Coupler now only supports the coupling of a scalar field or a field on a grid related to a horizontal grid that is decomposed in parallelization of a model. If you want to couple more kinds of fields, please contact us (liuli-cess@tsinghua.edu.cn)", field_name, annotation);

    if (buf_mark != BUF_MARK_IO_FIELD_REG)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, buf_mark >= 0, "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the parameter of the mark (\"buf_mark\") of the field instance cannot be a negative integer (currently is %d). Please check the model code with the annotation \"%s\"",
                         field_name, buf_mark, annotation);

    existing_field_instance_instance = search_field_instance(field_name, decomp_id, comp_or_grid_id, buf_mark);
    if (existing_field_instance_instance != NULL)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_field_instance/CCPL_start_chunk_field_instance_registration\" to register a field instance: cannot register an instance of the field of \"%s\" again (the corresponding annotation is \"%s\") because this field instance has been registered before (the corresponding annotation is \"%s\")", 
                         field_name, annotation, annotation_mgr->get_annotation(existing_field_instance_instance->get_field_instance_id(), "allocate field instance"));

    new_field_instance = new Field_mem_info(field_name, decomp_id, comp_or_grid_id, buf_mark, unit, data_type, annotation, (buf_mark!=BUF_MARK_IO_FIELD_REG) && (usage_tag&REG_FIELD_TAG_CPL) == REG_FIELD_TAG_CPL);
	if (data_type != NULL && new_field_instance->get_size_of_field() > 0)
	    EXECUTION_REPORT(REPORT_ERROR, comp_id, field_size == new_field_instance->get_size_of_field(), "Error happens when calling the API \"CCPL_register_field_instance\" to register a field instance of \"%s\": the size of the model data buffer (currently is %d) is different from the size determined by the parallel decomposition and grid (currently is %ld). Please check the model code with the annotation \"%s\"",
    	                 field_name, field_size, new_field_instance->get_size_of_field(), annotation);
    new_field_instance->set_field_instance_id(TYPE_FIELD_INST_ID_PREFIX|fields_mem.size(), annotation);
	new_field_instance->set_usage_tag(usage_tag);
	if (data_type != NULL)
	    new_field_instance->reset_mem_buf(data_buffer, true, usage_tag);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, usage_tag >= 0 && usage_tag <= 3, "Error happens when calling the API \"CCPL_register_field_instance/CCPL_start_chunk_field_instance_registration\" to register a field instance of \"%s\": the value of the parameter \"usage_tag\" (%d) is wrong. The right value should be between 1 and 3. Please check the model code with the annotation \"%s\"", field_name, usage_tag, annotation);
    fields_mem.push_back(new_field_instance);

    return new_field_instance->get_field_instance_id();
}


bool Memory_mgt::check_is_legal_field_instance_id(int field_instance_id)
{
    if ((field_instance_id&TYPE_ID_PREFIX_MASK) != TYPE_FIELD_INST_ID_PREFIX)
        return false;

    return (field_instance_id&TYPE_ID_SUFFIX_MASK) < fields_mem.size();
}


Field_mem_info *Memory_mgt::get_field_instance(int field_instance_id)
{
    if (!check_is_legal_field_instance_id(field_instance_id))
        return NULL;

    return fields_mem[field_instance_id&TYPE_ID_SUFFIX_MASK];
}


void Memory_mgt::copy_field_data_values(Field_mem_info *dst_field_inst, Field_mem_info *src_field_inst)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_field_inst->get_decomp_id() == src_field_inst->get_decomp_id() && dst_field_inst->get_comp_id() == src_field_inst->get_comp_id() && dst_field_inst->get_grid_id() == src_field_inst->get_grid_id() && words_are_the_same(src_field_inst->get_data_type(), dst_field_inst->get_data_type()),
                     "Software erorr in Memory_mgt::copy_field_data_values");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_field_inst->get_data_buf() != NULL && src_field_inst->get_data_buf() != NULL || dst_field_inst->get_data_buf() == NULL && src_field_inst->get_data_buf() == NULL, "Software error in Memory_mgt::copy_field_data_values");
	if (dst_field_inst->get_data_buf() != NULL)
	    memcpy(dst_field_inst->get_data_buf(), src_field_inst->get_data_buf(), dst_field_inst->get_size_of_field()*get_data_type_size(dst_field_inst->get_data_type()));

	for (int i = 0; i < dst_field_inst->get_num_chunks(); i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_field_inst->get_chunk_buf(i) != NULL && src_field_inst->get_chunk_buf(i) != NULL || dst_field_inst->get_chunk_buf(i) == NULL && src_field_inst->get_chunk_buf(i) == NULL, "Software error in Memory_mgt::copy_field_data_values");
	    memcpy(dst_field_inst->get_chunk_buf(i), src_field_inst->get_chunk_buf(i), dst_field_inst->get_chunk_data_buf_size(i)*get_data_type_size(dst_field_inst->get_data_type()));
	}

	dst_field_inst->define_field_values(false);
}


void Memory_mgt::get_comp_existing_registered_field_insts(std::vector<Field_mem_info *> &existing_field_insts, int comp_id)
{
	existing_field_insts.clear();
	for (int i = 0; i < fields_mem.size(); i ++)
		if (fields_mem[i]->get_comp_id() == comp_id && fields_mem[i]->get_is_registered_model_buf())
			existing_field_insts.push_back(fields_mem[i]);
}

