/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/

#include "global_data.h"
#include "runtime_datatype_transformer.h"
#include "cor_global_data.h"


Runtime_datatype_transformer::Runtime_datatype_transformer(Field_mem_info *src_field, Field_mem_info *dst_field)
{
    src_fields.push_back(src_field);
    dst_fields.push_back(dst_field);
}


bool Runtime_datatype_transformer::run(bool bypass_timer)
{
    transform_fields_datatype();
    
    return true;
}


void Runtime_datatype_transformer::handler_datatype_transformation_of_array(const char *data_type_dst, const char *data_type_src, void *data_dst, const void *data_src, const int num_local_cells)
{
	if (words_are_the_same(data_type_src, DATA_TYPE_DOUBLE)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_FLOAT))
			transform_datatype_of_arrays((double*)data_src, (float*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error1 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else if (words_are_the_same(data_type_src, DATA_TYPE_FLOAT)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_DOUBLE))
			transform_datatype_of_arrays((float*)data_src, (double*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error2 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else if (words_are_the_same(data_type_src, DATA_TYPE_LONG)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
			transform_datatype_of_arrays((long*)data_src, (int*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
			transform_datatype_of_arrays((long*)data_src, (short*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
			transform_datatype_of_arrays((long*)data_src, (bool*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error3 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else if (words_are_the_same(data_type_src, DATA_TYPE_INT)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
			transform_datatype_of_arrays((int*)data_src, (long*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
			transform_datatype_of_arrays((int*)data_src, (short*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
			transform_datatype_of_arrays((int*)data_src, (bool*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error4 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else if (words_are_the_same(data_type_src, DATA_TYPE_SHORT)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
			transform_datatype_of_arrays((short*)data_src, (long*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
			transform_datatype_of_arrays((short*)data_src, (int*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_BOOL))
			transform_datatype_of_arrays((short*)data_src, (bool*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error5 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else if (words_are_the_same(data_type_src, DATA_TYPE_BOOL)) {
		if (words_are_the_same(data_type_dst, DATA_TYPE_LONG))
			transform_datatype_of_arrays((bool*)data_src, (long*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_INT))
			transform_datatype_of_arrays((bool*)data_src, (int*)data_dst, num_local_cells);
		else if (words_are_the_same(data_type_dst, DATA_TYPE_SHORT))
			transform_datatype_of_arrays((bool*)data_src, (short*)data_dst, num_local_cells);
		else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error6 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");
	}
	else EXECUTION_REPORT(REPORT_ERROR,-1, "C-Coupler software error7 in Runtime_datatype_transformer::handler_datatype_transformation_of_array");

}


void Runtime_datatype_transformer::transform_fields_datatype()
{
    char *data_type_src, *data_type_dst;


    for (int i = 0; i < src_fields.size(); i ++) {
		src_fields[i]->check_field_sum(report_internal_log_enabled, true, "(src value) before data type transformation");
        src_fields[i]->use_field_values("");
        dst_fields[i]->define_field_values(false);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_fields[i]->get_num_chunks() == dst_fields[i]->get_num_chunks(), "Software error in Runtime_datatype_transformer::transform_fields_datatype");
        data_type_src = src_fields[i]->get_field_data()->get_grid_data_field()->data_type_in_application;
        data_type_dst = dst_fields[i]->get_field_data()->get_grid_data_field()->data_type_in_application;
		if (src_fields[i]->get_num_chunks() == 0)
			handler_datatype_transformation_of_array(data_type_dst, data_type_src, dst_fields[i]->get_data_buf(), src_fields[i]->get_data_buf(), src_fields[i]->get_field_data()->get_grid_data_field()->required_data_size);
		else for (int j = 0; j < src_fields[i]->get_num_chunks(); j ++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_fields[i]->get_chunk_data_buf_size(j) == src_fields[i]->get_chunk_data_buf_size(j), "Software error in Runtime_datatype_transformer::transform_fields_datatype");
			handler_datatype_transformation_of_array(data_type_dst, data_type_src, dst_fields[i]->get_chunk_buf(j), src_fields[i]->get_chunk_buf(j), dst_fields[i]->get_chunk_data_buf_size(j));
		}
		dst_fields[i]->check_field_sum(report_internal_log_enabled, true, "(dst value) after data type transformation");
    }
}

