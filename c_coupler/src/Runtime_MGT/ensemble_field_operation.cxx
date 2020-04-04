/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "ensemble_field_operation.h"


/* This implementation does not take consideration of mask */
template <class T> void member_to_set_operation_template(T **member_fields_data_buffer, T *set_field_data_buffer, int num_members, int field_size, int operation, int specified_member_index)
{
	if (operation == ENSEMBLE_OP_TYPE_ANY)
		memcpy(set_field_data_buffer, member_fields_data_buffer[specified_member_index], field_size*sizeof(T));
	else {
		for (int i = 0; i < field_size; i ++) {
			T ensemble_result = (T) 0;
			if (operation == ENSEMBLE_OP_TYPE_MAX || operation == ENSEMBLE_OP_TYPE_MIN)
				ensemble_result = member_fields_data_buffer[0][i];
			for (int j = 0; j < num_members; j ++) {
				if (operation == ENSEMBLE_OP_TYPE_SUM || operation == ENSEMBLE_OP_TYPE_MEAN || operation == ENSEMBLE_OP_TYPE_ANOMALY)
					ensemble_result += member_fields_data_buffer[j][i];
				else if (operation == ENSEMBLE_OP_TYPE_MIN) {
					if (ensemble_result > member_fields_data_buffer[j][i])
						ensemble_result = member_fields_data_buffer[j][i];
				}
				else if (operation == ENSEMBLE_OP_TYPE_MAX) {
					if (ensemble_result < member_fields_data_buffer[j][i])
						ensemble_result = member_fields_data_buffer[j][i];					
				}
				else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in member_to_set_operation_template");
			}
			if (operation == ENSEMBLE_OP_TYPE_MEAN)
				set_field_data_buffer[i] = ensemble_result / (T)num_members;
			else if (operation == ENSEMBLE_OP_TYPE_ANOMALY)
				set_field_data_buffer[i] = set_field_data_buffer[i] - ensemble_result/(T)num_members;
			else set_field_data_buffer[i] = ensemble_result;
		}
	}
}


Member_to_set_operation::Member_to_set_operation(std::vector<Field_mem_info*> &member_fields_inst, Field_mem_info *set_field_inst, int operation, int specified_member_index)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, member_fields_inst.size() > 0 && operation >= ENSEMBLE_OP_TYPE_SUM && operation <= ENSEMBLE_OP_TYPE_ANY, "Software error in Member_to_set_operation::Member_to_set_operation");
	if (operation == ENSEMBLE_OP_TYPE_ANY)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, specified_member_index >= 0 && specified_member_index < member_fields_inst.size(), "Software error in Member_to_set_operation::Member_to_set_operation");
	this->specified_member_index = specified_member_index;
	for (int i = 0; i < member_fields_inst.size(); i ++)
		this->member_fields_inst.push_back(member_fields_inst[i]);
	this->set_field_inst = set_field_inst;
	member_fields_data_buffer = new void *[member_fields_inst.size()];
	if (set_field_inst != NULL)
		set_field_data_buffer = set_field_inst->get_data_buf();
	else set_field_data_buffer = NULL;
	operation_type = operation;
	field_size = member_fields_inst[0]->get_size_of_field();
	for (int i = 0; i < member_fields_inst.size(); i ++) {
		if (set_field_inst != NULL)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, member_fields_inst[i]->get_data_buf() != set_field_inst->get_data_buf() && member_fields_inst[i]->get_comp_id() == set_field_inst->get_comp_id() && member_fields_inst[i]->get_grid_id() == set_field_inst->get_grid_id() && member_fields_inst[i]->get_decomp_id() == set_field_inst->get_decomp_id() && words_are_the_same(member_fields_inst[i]->get_data_type(), set_field_inst->get_data_type()));
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, member_fields_inst[i]->get_comp_id() == member_fields_inst[0]->get_comp_id() && member_fields_inst[i]->get_grid_id() == member_fields_inst[0]->get_grid_id() && member_fields_inst[i]->get_decomp_id() == member_fields_inst[0]->get_decomp_id() && words_are_the_same(member_fields_inst[i]->get_data_type(), member_fields_inst[0]->get_data_type()));
		member_fields_data_buffer[i] = member_fields_inst[i]->get_data_buf();
	}
//	if (operation != ENSEMBLE_OP_TYPE_ANY)
//		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(set_field_inst->get_data_type(), DATA_TYPE_DOUBLE) || words_are_the_same(set_field_inst->get_data_type(), DATA_TYPE_FLOAT), "Software error in Member_to_set_operation::Member_to_set_operation");
	data_type = strdup(member_fields_inst[0]->get_data_type());
}


void Member_to_set_operation::execute()
{
	if (words_are_the_same(data_type, DATA_TYPE_BOOL))
		member_to_set_operation_template((bool**)member_fields_data_buffer, (bool*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_CHAR))
		member_to_set_operation_template((char**)member_fields_data_buffer, (char*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
		member_to_set_operation_template((double**)member_fields_data_buffer, (double*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
		member_to_set_operation_template((float**)member_fields_data_buffer, (float*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_INT))
		member_to_set_operation_template((int**)member_fields_data_buffer, (int*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_LONG))
		member_to_set_operation_template((long**)member_fields_data_buffer, (long*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else if (words_are_the_same(data_type, DATA_TYPE_SHORT))
		member_to_set_operation_template((short**)member_fields_data_buffer, (short*) set_field_data_buffer, member_fields_inst.size(), field_size, operation_type, specified_member_index);
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Member_to_set_operation::execute");
	set_field_inst->check_field_sum(report_internal_log_enabled, true, "Set field after member_to_set_operation");
}


Member_to_set_operation::~Member_to_set_operation()
{
	delete [] member_fields_data_buffer;
	delete [] data_type;
}


