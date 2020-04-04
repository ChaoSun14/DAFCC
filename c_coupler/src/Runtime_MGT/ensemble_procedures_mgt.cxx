/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include <dlfcn.h>
#include "global_data.h"
#include "ensemble_procedures_mgt.h"
#include "ensemble_field_operation.h"
#include <algorithm>
#include <unistd.h>
#include <cstring>

Field_instances_operation::Field_instances_operation(Ensemble_procedures_inst *ensemble_procedures_inst, TiXmlElement *field_instance_XML_element)
{
	int line_number;
	int member_id;

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start Field_instances_operation::Field_instances_operation");
	this->ensemble_procedures_inst = ensemble_procedures_inst;
	if (is_XML_setting_on(ensemble_procedures_inst->get_member_comp_id(), field_instance_XML_element, ensemble_procedures_inst->get_XML_file_name(), "the status \"field_instances\"", "ensemble procedures configuration")) {
		this->do_field_instances_operation = true;
		this->use_statistical_method = false;
		this->do_ensemble_op = false;
        const char *field_instances_statistical_method = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_instance_XML_element, "statistical_method", ensemble_procedures_inst->get_XML_file_name(), line_number, "default statistical method of all field instances", "ensemble procedures configuration", true);
        	check_statistical_method_format(field_instances_statistical_method, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of default statistical method of all field instances");
        const char *field_instances_ensemble_op = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_instance_XML_element, "ensemble_op", ensemble_procedures_inst->get_XML_file_name(), line_number, "default ensemble operation of all field instances", "ensemble procedures configuration", true);   
        	this->field_instances_member_id = check_ensemble_op_format(field_instances_ensemble_op, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of default ensemble operation of all field instances");
        this->field_instances_statistical_method = strdup(field_instances_statistical_method);
        this->field_instances_ensemble_op = strdup(field_instances_ensemble_op);
        for (TiXmlNode *field_XML_node = field_instance_XML_element->FirstChildElement(); field_XML_node != NULL; field_XML_node = field_XML_node->NextSibling()) {
        	if (field_XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        	TiXmlElement *field_XML_element = field_XML_node->ToElement();
        	EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(field_XML_element->Value(),"field"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"field\". Please verify the XML file arround the line number %d.", ensemble_procedures_inst->get_XML_file_name(), field_instance_XML_element->Row());
        	const char *field_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_XML_element, "name", ensemble_procedures_inst->get_XML_file_name(), line_number, "name of a field", "configuration of the attributes of shared fields for ensemble procedures", true);
        	bool field_existed = true;
        	for (int i = 0; i < ensemble_procedures_inst->get_num_specified_field_instances(); i ++) {
        		//EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "API_specified_field_insts[%d]->get_field_name(): %s", i+1, ensemble_procedures_inst->API_specified_field_insts[i]->get_field_name());
        		if (words_are_the_same(ensemble_procedures_inst->API_specified_field_insts[i]->get_field_name(), field_name)) {
        			field_existed = true;
        			break;
        		}
        		else field_existed = false;
        	}
            EXECUTION_REPORT(REPORT_ERROR, -1, field_existed, "The field name \"%s\" in the XML configuration file \"%s\" has not been registered in C-Coupler. Please verify the XML file arround the line number %d.", field_name, ensemble_procedures_inst->get_XML_file_name(), field_XML_element->Row());
        	const field_op *existing_field = search_field_op(field_name);
        	if (existing_field != NULL) {
        		EXECUTION_REPORT(REPORT_ERROR, -1, false, "Cannot spefify the attributes of field \"%s\" in the XML file \"%s\" around the line number %d again because it has already been specified around the line number %d", field_name, ensemble_procedures_inst->get_XML_file_name(), line_number, existing_field->line_number);
        	}
        	const char *field_statistical_method = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_XML_element, "statistical_method", ensemble_procedures_inst->get_XML_file_name(), line_number, "statistical method of a field", "configuration of the attributes of shared fields for ensemble procedures", false);
        	if (field_statistical_method != NULL) {
        		check_statistical_method_format(field_statistical_method, field_name, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of statistical method of the field instance");
        	}
        	const char *field_ensemble_op = get_XML_attribute(-1, CCPL_NAME_STR_LEN, field_XML_element, "ensemble_op", ensemble_procedures_inst->get_XML_file_name(), line_number, "ensemble operation of a field", "configuration of the attributes of shared fields for ensemble procedures", false);
        	if (field_ensemble_op != NULL) {
        		member_id = check_ensemble_op_format(field_ensemble_op, field_name, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of ensemble operation of the field instance");
        	}
        	add_field_op(field_name, field_statistical_method, field_ensemble_op, member_id, line_number);
        }  
    }
    else {
      	this->do_field_instances_operation = false;
      	this->use_statistical_method = false;
      	this->do_ensemble_op = false;
       	EXECUTION_REPORT(REPORT_LOG, ensemble_procedures_inst->get_member_comp_id(), true, "The status of \"field_instances\" in the XML configuration file \"%s\" is \"off\", which means the field instances operation (\"statistical_method\" and \"ensemble_op\") is not used in the ensemble procedures instance.", ensemble_procedures_inst->get_XML_file_name());
    }
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish Field_instances_operation::Field_instances_operation");
}

const field_op *Field_instances_operation::search_field_op(const char *field_name)
{
    for (int i = 0; i < fields_op.size(); i ++)
        if (words_are_the_same(field_name, fields_op[i].field_name))
            return &fields_op[i];
    return NULL;
}

void Field_instances_operation::add_field_op(const char *field_name, const char *field_statistical_method, const char *field_ensemble_op, int member_id, int line_number)
{
    field_op local_op;

    strcpy(local_op.field_name, field_name);
    strcpy(local_op.field_statistical_method, field_statistical_method);
    strcpy(local_op.field_ensemble_op, field_ensemble_op);
    local_op.member_id = member_id;
    local_op.line_number = line_number;
    fields_op.push_back(local_op);
    EXECUTION_REPORT(REPORT_ERROR, -1, search_field_op(local_op.field_name) == &(fields_op[fields_op.size()-1]), "Software error in Field_instances_operation::add_field_op");
}

void Field_instances_operation::check_statistical_method_format(const char *statistical_method, const char *XML_file_name, int line_number, const char *statistical_method_annotation)
{
	if (words_are_the_same(statistical_method, STATISTICAL_METHOD_INST)) this->use_statistical_method = false;
	else{
		if (words_are_the_same(statistical_method, STATISTICAL_METHOD_AVER) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MAX) || words_are_the_same(statistical_method, STATISTICAL_METHOD_ACCUM) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MIN)) this->use_statistical_method = true;
		else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of statistical method \"%s\" with the annotation \"%s\" is wrong (must be one of \"inst/aver/accum/max/min\"). Please verify the XML file \"%s\" arround the line number %d.", statistical_method, statistical_method_annotation, XML_file_name, line_number);
	}
}

void Field_instances_operation::check_statistical_method_format(const char *statistical_method, const char *field_name, const char *XML_file_name, int line_number, const char *statistical_method_annotation)
{
	if (words_are_the_same(statistical_method, STATISTICAL_METHOD_INST)) this->use_statistical_method = false;
	else{
		if (words_are_the_same(statistical_method, STATISTICAL_METHOD_AVER) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MAX) || words_are_the_same(statistical_method, STATISTICAL_METHOD_ACCUM) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MIN)) this->use_statistical_method = true;
		else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of statistical method \"%s\" for the field \"%s\" with the annotation \"%s\" is wrong (must be one of \"inst/aver/accum/max/min\"). Please verify the XML file \"%s\" arround the line number %d.", statistical_method, field_name, statistical_method_annotation, XML_file_name, line_number);
	}
}

int Field_instances_operation::check_ensemble_op_format(const char *ensemble_op, const char *XML_file_name, int line_number, const char *ensemble_op_annotation)
{
	int member_id;
    //char member_str[5]; 
	member_id = -1;
	if (strstr(ensemble_op, "mem_") == NULL) {
		if (words_are_the_same(ensemble_op, ENSEMBLE_OP_NONE)) this->do_ensemble_op = false;
		else{
			if (words_are_the_same(ensemble_op, ENSEMBLE_OP_GATHER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_ANOM)) {
                this->do_ensemble_op = true;
                this->use_set_grids = true;
            } 
            else {
                if (words_are_the_same(ensemble_op, ENSEMBLE_OP_AVER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MAX) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MIN)) {
                    this->do_ensemble_op = true;
                    this->use_set_grids = false;
                }
                else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" with the annotation \"%s\" is wrong (must be one of \"none/gather/aver/anom/max/min/mem_%d\"). Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, ensemble_op_annotation, XML_file_name, line_number);
            }
		}
	}
	else{
        if (sscanf(ensemble_op, "mem_%d", &member_id) != 1) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" with the annotation \"%s\" is wrong (must be one of \"none/gather/aver/anom/max/min/mem_%d\"). Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, ensemble_op_annotation, XML_file_name, line_number);
		if (member_id == 0 || member_id > ensemble_procedures_inst->num_ens_members) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" with the annotation \"%s\" is wrong, the ensemble member ID should >= 1 and <= %d. Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, ensemble_op_annotation, ensemble_procedures_inst->num_ens_members, XML_file_name, line_number); 
		this->do_ensemble_op = true;
        this->use_set_grids = false;
	}
	return member_id;
}

int Field_instances_operation::check_ensemble_op_format(const char *ensemble_op, const char *field_name, const char *XML_file_name, int line_number, const char *ensemble_op_annotation)
{
	int member_id;
	member_id = -1;
	if (do_ensemble_op){
		if (words_are_the_same(ensemble_op, ENSEMBLE_OP_NONE)) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" for the field \"%s\" with the annotation \"%s\" conflicts with default values \"%s\". Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, field_name, ensemble_op_annotation, this->field_instances_ensemble_op, XML_file_name, line_number);
		if (strstr(ensemble_op, "mem_") == NULL) {
			if (words_are_the_same(ensemble_op, ENSEMBLE_OP_GATHER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_ANOM)) {
                this->do_ensemble_op = true;
                this->use_set_grids = true;
            } 
            else {
                if (words_are_the_same(ensemble_op, ENSEMBLE_OP_AVER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MAX) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MIN)) {
                    this->do_ensemble_op = true;
                    this->use_set_grids = false;
                }
                else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" with the annotation \"%s\" is wrong (must be one of \"none/gather/aver/anom/max/min/mem_%d\"). Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, ensemble_op_annotation, XML_file_name, line_number);
            }
		}
		else{
			if (sscanf(ensemble_op, "mem_%d", &member_id) == 0) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" for the field \"%s\" with the annotation \"%s\" is wrong (must be one of \"none/gather/aver/anom/max/min/mem_%d\"). Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, field_name, ensemble_op_annotation, XML_file_name, line_number);
			if (member_id == 0 || member_id > ensemble_procedures_inst->num_ens_members) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" for the field \"%s\" with the annotation \"%s\" is wrong, the ensemble member ID should >= 1 and <= %d. Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, field_name, ensemble_op_annotation, ensemble_procedures_inst->num_ens_members, XML_file_name, line_number); 
			this->do_ensemble_op = true;
            this->use_set_grids = false;
		}
	}
	else if(!words_are_the_same(ensemble_op, ENSEMBLE_OP_NONE)) EXECUTION_REPORT(REPORT_ERROR, -1, false, "The attributes of ensemble operation \"%s\" for the field \"%s\" with the annotation \"%s\" conflicts with default values \"%s\". Please verify the XML file \"%s\" arround the line number %d.", ensemble_op, field_name, ensemble_op_annotation, this->field_instances_ensemble_op, XML_file_name, line_number);	
	return member_id;	
}
const char *Field_instances_operation::get_field_op_statistical_method(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_statistical_method;
    
    return search_field_op(field_name)->field_statistical_method;
}
const char *Field_instances_operation::get_field_op_ensemble_op(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_ensemble_op;
    
    return search_field_op(field_name)->field_ensemble_op;
}
int Field_instances_operation::get_field_op_member_id(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_member_id;
    
    return search_field_op(field_name)->member_id;
}

Ensemble_procedures_inst::Ensemble_procedures_inst(int instance_id, const char *inst_name, int set_comp_id, int member_comp_id, int size_field_inst, int size_grids, int size_decomps, int size_timers, int size_controls,  
                                                   const int *field_inst_ids, const int *grid_ids, const int *decomp_ids, const int *timer_ids, const int *control_vars, const char *annotation)
{
	char tmp_string[NAME_STR_SIZE];
	int line_number;
	int tmp_int, num_total_active_procs;

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start Ensemble_procedures_inst");
	check_and_verify_name_format_of_string_for_API(member_comp_id, inst_name, API_ID_ENSEMBLE_PROC_INST_INIT, "the name of the ensemble procedure(s) instance", annotation);
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(member_comp_id, "Ensemble_procedures_inst::Ensemble_procedures_inst");
    synchronize_comp_processes_for_API(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", annotation);
    check_API_parameter_string(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", inst_name, "inst_name", annotation);
    
	check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_controls, "array size of \"control_vars\" (array size)", annotation);	
	check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_field_inst, "array size of \"field_inst_ids\" (array size)", annotation);
	check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_decomps, "array size of \"decomp_ids\" (array size)", annotation);
	check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_grids, "array size of \"grid_ids\" (array size)", annotation);
	check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_timers, "array size of \"timer_ids\" (array size)", annotation);
	check_API_parameter_data_array(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_controls, sizeof(int), (const char*) control_vars, "control_vars", annotation);
	for (int i = 0; i < size_controls; i ++)
		this->control_vars.push_back(control_vars[i]);
    for (int i = 0; i < size_field_inst; i ++) {
        check_API_parameter_field_instance(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", field_inst_ids[i], "field_inst_ids", annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, memory_manager->get_field_instance(field_inst_ids[i])->get_comp_id() == member_comp_id, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d field instance specified by the input parameter \"field_ids\"\" does not belong to the component model corresponding to the input parameter of \"member_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
		this->API_specified_field_insts.push_back(memory_manager->get_field_instance(field_inst_ids[i]));
    }
	for (int i = 0; i < size_timers; i ++) {
        check_API_parameter_timer(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", timer_ids[i], "timer_ids", annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, timer_mgr->get_timer(timer_ids[i])->get_comp_id() == member_comp_id, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d timer specified by the input parameter \"timer_ids\"\" does not belong to the component model corresponding to the input parameter of \"member_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
		this->timer_ids.push_back(timer_ids[i]);
	}	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, size_grids == size_decomps, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the array size of the input parameter \"grid_ids\" (%d) is different from the array size of the input parameter \"decomp_ids\" (%d). Please verify at the model code with the annotation \"%s\".", inst_name, size_grids, size_decomps, annotation);
	for (int i = 0; i < size_grids; i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, original_grid_mgr->is_grid_id_legal(grid_ids[i]), "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d element of the input parameter \"grid_ids\" is not a legal ID of grids Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, original_grid_mgr->search_grid_info(grid_ids[i])->get_comp_id() == member_comp_id, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d grid specified by the input parameter \"grid_ids\"\" does not belong to the component model corresponding to the input parameter of \"member_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);				
		check_API_parameter_string(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", original_grid_mgr->search_grid_info(grid_ids[i])->get_grid_name(), "\"grid_ids\"", annotation);
		if (decomp_ids[i] == -1 || decomp_ids[i] == CCPL_NULL_INT)
			tmp_int = -1;
		else tmp_int = 0;
		check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", tmp_int, "\"decomp_ids\"", annotation);	
		if (decomp_ids[i] != -1 && decomp_ids[i] != CCPL_NULL_INT) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, decomps_info_mgr->is_decomp_id_legal(decomp_ids[i]), "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d element of the input parameter \"decomp_ids\" is not a legal ID of parallel decomposition. Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_comp_id() == member_comp_id, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d parallel decomposition specified by the input parameter \"decomp_ids\"\" does not belong to the component model corresponding to the input parameter of \"member_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, original_grid_mgr->search_grid_info(decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_grid_id())->get_original_CoR_grid()->is_subset_of_grid(original_grid_mgr->search_grid_info(grid_ids[i])->get_original_CoR_grid()), "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d parallel decomposition in the input parameter \"decomp_ids\"\" is not consistent with the corresponding grid speicified in the input parameter \"grid_ids\". Please verify at the model code with the annotation \"%s\".", inst_name, i+1, annotation);
			check_API_parameter_string(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_decomp_name(), "\"decomp_ids\"", annotation);
		}
		this->grid_ids.push_back(grid_ids[i]);
		this->decomp_ids.push_back(decomp_ids[i]);
	}
    EXECUTION_REPORT(REPORT_ERROR, member_comp_id, components_time_mgrs->get_time_mgr(set_comp_id)->get_time_step_in_second() > 0, "Error happers when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the time step of the corresponding component model \"%s\" has not been set yet. Please specify the time step before the model code with the annotation \"%s\"", inst_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(set_comp_id, true, "Ensemble_procedures_inst::Ensemble_procedures_inst")->get_comp_full_name(), annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish API parameters check in Ensemble_procedures_inst");
	this->set_comp_id = set_comp_id;
	this->member_comp_id = member_comp_id;
	this->instance_name = strdup(inst_name);
	this->instance_id = instance_id;
	this->local_comm = comm;
    int member_id;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true,"%s", comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name());
    if (strstr(comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name(),"_member") != NULL) { 
        if (sscanf(strstr(comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name(),"_member"), "_member%d", &member_id) !=1 ) 
        EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to get the ensemble id of current process");
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to get the ensemble id of current process, the componant name should be extented with \"_member[id]\", while the current name is \"%s\"",comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name());
    this->proc_member_id = member_id;
    this->num_ens_members = comp_comm_group_mgt_mgr->get_num_members_in_ensemble_set(comp_comm_group_mgt_mgr->search_global_node(this->set_comp_id), comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id));
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The number of ensemble members is %d, and the ensemble id of current process is %d", this->num_ens_members, this->proc_member_id);

	MPI_Comm_rank(this->local_comm, &proc_id_in_local_comm);
	MPI_Comm_size(this->local_comm, &num_proc_in_local_comm);
	
	annotation_mgr->add_annotation(instance_id, "registering an instance of ensemble procedures", annotation);

	for (int i = 0; i < API_specified_field_insts.size(); i ++)
		for (int j = i+1; j < API_specified_field_insts.size(); j ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, !words_are_the_same(API_specified_field_insts[i]->get_field_name(),API_specified_field_insts[j]->get_field_name()), "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the input parameter \"field_inst_ids\" contains multiple elements corresponding to the same field \"%s\" (at number %d and %d elements). Please verify at the model code with the annotation \"%s\" (please specify at most one element for this field).", inst_name, API_specified_field_insts[i]->get_field_name(), i+1, j+1, annotation);
	
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations in Ensemble_procedures_inst");
	this->XML_file_name = NULL;
	this->procedures_name = NULL;
	this->procedures_type = NULL;
	this->procedures_dl_name = NULL;
	TiXmlElement *config_XML_element = get_XML_file_with_configuration();
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, config_XML_element != NULL, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\" of ensemble procedures package at the model code with the annotation \"%s\": no corresponding configuration information can be found from the corresponding XML files under the directory \"%s\". Please verify.", inst_name, annotation, comp_comm_group_mgt_mgr->get_ensemble_procedure_config_dir());	
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about external_procedures in Ensemble_procedures_inst");
	//TiXmlNode *field_XML_node = config_XML_element->FirstChildElement(); 
    TiXmlElement *external_procedures_XML_element = config_XML_element->FirstChildElement();
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(external_procedures_XML_element->Value(),"external_procedures"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"external_procedures\". Please verify the XML file arround the line number %d.", XML_file_name, external_procedures_XML_element->Row());
    if (is_XML_setting_on(member_comp_id, external_procedures_XML_element, XML_file_name, "the status \"external_procedures\"", "ensemble procedures configuration")) {
      	const char *p_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "procedures_name", XML_file_name, line_number, "name of external procedures", "ensemble procedures configuration", true);
       	const char *p_type = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "type", XML_file_name, line_number, "type of external procedures", "ensemble procedures configuration", true);
       	const char *p_dl_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "dl_name", XML_file_name, line_number, "dynamic library name of external procedures", "ensemble procedures configuration", true);
       	//EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "procedures_name, procedures_type, procedures_dl_name, %s, %s, %s", p_name, p_type, p_dl_name);
       	this->procedures_name = strdup(p_name);
       	this->procedures_type = strdup(p_type);
       	this->procedures_dl_name = strdup(p_dl_name);
    } 
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "The status of \"external_procedures\" in the XML configuration file \"%s\" should be \"on\". Please verify the XML file arround the line number %d.", XML_file_name, external_procedures_XML_element->Row());
    
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about periodic_timer in Ensemble_procedures_inst");
    TiXmlElement *periodic_timer_XML_element = external_procedures_XML_element->NextSiblingElement();
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(periodic_timer_XML_element->Value(),"periodic_timer"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"periodic_timer\". Please verify the XML file arround the line number %d.", XML_file_name, periodic_timer_XML_element->Row());
    if (is_XML_setting_on(member_comp_id, periodic_timer_XML_element, XML_file_name, "the status \"periodic_timer\"", "ensemble procedures configuration")) {
    	const char *period_unit = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "period_unit", XML_file_name, line_number, "period unit of periodic timer", "ensemble procedures configuration", false);
       	if (period_unit != NULL) {
        	const char *period_count = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "period_count", XML_file_name, line_number, "period count of periodic timer", "ensemble procedures configuration", true);
        	const char *local_lag_count = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "local_lag_count", XML_file_name, line_number, "local lag count of periodic timer", "ensemble procedures configuration", true);
       		this->use_periodic_timer = true;
        	this->periodic_timer_id = timer_mgr->define_timer(set_comp_id, period_unit, atoi(period_count), atoi(local_lag_count), 0, "define periodic timer of ensemble procedures instance");
        }
        else this->use_periodic_timer = false;     	
    } 
    else {
        this->use_periodic_timer = false;
        EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"periodic_timer\" in the XML configuration file \"%s\" is \"off\", which means the periodic timer is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about field_instance in Ensemble_procedures_inst");
    TiXmlElement *field_instance_XML_element = periodic_timer_XML_element->NextSiblingElement();
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(field_instance_XML_element->Value(),"field_instances"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"field_instances\". Please verify the XML file arround the line number %d.", XML_file_name, field_instance_XML_element->Row());
    this->field_instances_op = new Field_instances_operation(this, field_instance_XML_element);
    
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about config_scripts in Ensemble_procedures_inst");
    TiXmlElement *config_scripts_XML_element = field_instance_XML_element->NextSiblingElement();
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(config_scripts_XML_element->Value(),"config_scripts"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"config_scripts\". Please verify the XML file arround the line number %d.", XML_file_name, config_scripts_XML_element->Row());
    if (is_XML_setting_on(member_comp_id, config_scripts_XML_element, XML_file_name, "the status \"config_scripts\"", "ensemble procedures configuration")) {
        TiXmlNode *before_instance_script_XML_node = config_scripts_XML_element->FirstChildElement();
        TiXmlElement *before_instance_script_XML_element = before_instance_script_XML_node->ToElement();
        if (is_XML_setting_on(member_comp_id, before_instance_script_XML_element, XML_file_name, "the status \"before_instance_script\"", "ensemble procedures configuration")) {
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(before_instance_script_XML_element->Value(),"before_instance_script"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"before_instance_script\". Please verify the XML file arround the line number %d.", XML_file_name, before_instance_script_XML_element->Row());
            const char *before_instance_script_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, before_instance_script_XML_element, "name", XML_file_name, line_number, "name of configuration script before instance run", "configuration of the attributes of shared fields for ensemble procedures", true);
            this->before_instance_script = strdup(before_instance_script_name);
        }
        else {
            this->before_instance_script = NULL;
            EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"before_instance_script\" in the XML configuration file \"%s\" is \"off\", which means the configuration script before instance run is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
        }
        TiXmlNode *after_instance_script_XML_node = before_instance_script_XML_node->NextSibling();
        TiXmlElement *after_instance_script_XML_element = after_instance_script_XML_node->ToElement();
        if (is_XML_setting_on(member_comp_id, after_instance_script_XML_element, XML_file_name, "the status \"after_instance_script\"", "ensemble procedures configuration")) {
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(after_instance_script_XML_element->Value(),"after_instance_script"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"after_instance_script\". Please verify the XML file arround the line number %d.", XML_file_name, after_instance_script_XML_element->Row());
            const char *after_instance_script_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, after_instance_script_XML_element, "name", XML_file_name, line_number, "name of configuration script after instance run", "configuration of the attributes of shared fields for ensemble procedures", true);
            this->after_instance_script = strdup(after_instance_script_name);
        }
        else {
            this->after_instance_script = NULL;
            EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"after_instance_script\" in the XML configuration file \"%s\" is \"off\", which means the configuration script after instance run is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
        }
    }
    else {
        EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"config_scripts\" in the XML configuration file \"%s\" is \"off\", which means the configuration scripts are not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about namelist in Ensemble_procedures_inst");
    TiXmlElement *namelist_XML_element = config_scripts_XML_element->NextSiblingElement();
    EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(namelist_XML_element->Value(),"namelist"), "The XML element for specifying the attributes of a public field in the XML configuration file \"%s\" should be named \"namelist\". Please verify the XML file arround the line number %d.", XML_file_name, namelist_XML_element->Row());
    ensemble_procedures_mgr->add_ensemble_procedures_inst(this);
    
    if (this->get_field_instances_op()->if_do_field_instances_operation() && !this->get_field_instances_op()->if_do_ensemble_op()) do_none_ensemble_op_initialize();
    if (this->get_field_instances_op()->if_do_field_instances_operation() && this->get_field_instances_op()->if_do_ensemble_op()) do_ensemble_op_initialize();

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish Ensemble_procedures_inst");
}



Ensemble_procedures_inst::~Ensemble_procedures_inst()
{
}


TiXmlElement *Ensemble_procedures_inst::get_XML_file_with_configuration()
{
    char current_XML_file_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *member_comp_node = comp_comm_group_mgt_mgr->search_global_node(member_comp_id);
    TiXmlDocument *current_XML_file;

    Comp_comm_group_mgt_node *current_comp_node = member_comp_node;
    if (words_are_the_same(current_comp_node->get_comp_type(), COMP_TYPE_PSEUDO_COUPLED)) return NULL;
    sprintf(current_XML_file_name, "%s/%s_DA_config.xml", comp_comm_group_mgt_mgr->get_ensemble_procedure_config_dir(), comp_comm_group_mgt_mgr->search_global_node(set_comp_id)->get_comp_name());
    current_XML_file = open_XML_file_to_read(member_comp_id, current_XML_file_name, member_comp_node->get_comm_group(), false);
    if (current_XML_file == NULL) return NULL;
    TiXmlElement *XML_element = current_XML_file->FirstChildElement();
    for (TiXmlNode *XML_node = get_XML_first_child_of_unique_root(member_comp_id, XML_file_name, current_XML_file); XML_node != NULL; XML_node = XML_node->NextSibling()) {
        if (XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        XML_element = XML_node->ToElement();
        if (!words_are_the_same(XML_element->Value(),"da_instance"))
            continue;
        if (!is_XML_setting_on(member_comp_id, XML_element, XML_file_name, "the status \"da_instance\"", "ensemble data assimilation procedures configuration"))
            continue;           
        const char *current_procedures_name = get_XML_attribute(member_comp_id, -1, XML_element, "name", XML_file_name, line_number, "the name of the ensemble procedures package", "ensemble data assimilation procedures configuration", true);
        if (words_are_the_same(current_procedures_name, this->instance_name)) {
            XML_file_name = strdup(current_XML_file_name);
            XML_file = current_XML_file;
            EXECUTION_REPORT_LOG(REPORT_LOG, member_comp_id, true, "In the process of calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\" of ensemble data assimilation procedures package: find the XML file \"%s\" with the configuration information of the ensemble procedures package", this->instance_name, XML_file_name);
            return XML_element;
        }
    }
    delete current_XML_file;

    return NULL;
}
void Ensemble_procedures_inst::do_ensemble_op(std::vector<Field_mem_info*> &field_insts_list, std::vector<Field_mem_info*> &set_field_insts)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to ensemble_procedures_inst::do_ensemble_op");
    //ensemble_to_set_operation = new Member_to_set_operation();
    int ensemble_operation, mem_id;
    char op_mem[NAME_STR_SIZE];
    std::vector<Field_mem_info*> tmp_member_field_insts;
    for (int i = 0; i < field_insts_list.size(); i ++) {
        if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_GATHER)){
            //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "do_ensemble_op: ENSEMBLE_OP_GATHER");
            set_field_insts.push_back(this->ensemble_set_field_insts_before_ens_op[i]);
            set_field_insts.back()->define_field_values(false);
        } 
        //else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_ANOM)) 
        //ensemble_operation=ENSEMBLE_OP_TYPE_ANOMALY;
        else{
            mem_id = this->field_instances_op->get_field_op_member_id(field_insts_list[i]->get_field_name()); 
            sprintf(op_mem, "mem_%d", mem_id);
            if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_AVER)) 
                ensemble_operation = ENSEMBLE_OP_TYPE_MEAN;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_MAX)) 
                ensemble_operation = ENSEMBLE_OP_TYPE_MAX;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_MIN)) 
                ensemble_operation = ENSEMBLE_OP_TYPE_MIN;
            //else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_SUM)) 
            //ensemble_operation = ENSEMBLE_OP_TYPE_SUM;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), op_mem)){
                //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "do_ensemble_op: ENSEMBLE_OP_MEM_%d", mem_id);
                //set_field_insts.push_back(memory_manager->get_field_instance(ensemble_member_field_insts_id[mem_id-1][i]));
                ensemble_operation = ENSEMBLE_OP_TYPE_ANY;
            }            
            else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Ensemble_procedures_inst::do_ensemble_op");
            tmp_member_field_insts.clear();
            for (int j = 0; j < this->ensemble_set_field_insts_one_member.size(); j ++){
                if (words_are_the_same(field_insts_list[i]->get_field_name(), this->ensemble_set_field_insts_one_member[j]->get_field_name())){
                    for (int n = 0; n < this->num_ens_members; n ++) 
                        tmp_member_field_insts.push_back(memory_manager->get_field_instance(ensemble_member_field_insts_id[n][j]));
                    Member_to_set_operation *ensemble_to_set_operation = new Member_to_set_operation(tmp_member_field_insts, this->ensemble_set_field_insts_one_member[j], ensemble_operation, mem_id);
                    ensemble_to_set_operation->execute();
                    set_field_insts.push_back(this->ensemble_set_field_insts_one_member[j]);
                    delete ensemble_to_set_operation;
                }
            }
        }
    }
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish ensemble_procedures_inst::do_ensemble_op");
}


void Ensemble_procedures_inst::do_none_ensemble_op_initialize()
{
	int *fields_id = new int [API_specified_field_insts.size()];
	for (int i = 0; i < API_specified_field_insts.size(); i ++) {
        const char *data_type = API_specified_field_insts[i]->get_data_type();
        mirror_API_specified_field_insts.push_back(memory_manager->alloc_mem(API_specified_field_insts[i], BUF_MARK_ENS_DATA_TRANSFER, this->instance_id, data_type, false));
        memory_manager->copy_field_data_values(mirror_API_specified_field_insts[i], API_specified_field_insts[i]);
        mirror_API_specified_field_insts[i]->define_field_values(false);
        fields_id[i] = mirror_API_specified_field_insts[i]->get_field_instance_id();
        //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "mirror_API_specified_field_insts[%d], %s, %d, %g", i, mirror_API_specified_field_insts[i]->get_field_name(), mirror_API_specified_field_insts[i]->get_size_of_field(), mirror_API_specified_field_insts[i]->get_data_buf());
        //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "fields_id[%d]=%d, tmp_fields_id[%d]=%d", i, API_specified_field_insts[i]->get_field_instance_id(), i, fields_id[i]);
    }
    int member_periodic_timer_id = timer_mgr->define_timer(this->member_comp_id, timer_mgr->get_timer(this->periodic_timer_id)->get_frequency_unit(), timer_mgr->get_timer(this->periodic_timer_id)->get_frequency_count(), timer_mgr->get_timer(this->periodic_timer_id)->get_local_lag_count(), timer_mgr->get_timer(this->periodic_timer_id)->get_remote_lag_count(), "define the periodic timer of ensemble procedures instance at do_none_ensemble_op_initialize");
    this->periodic_timer_id = member_periodic_timer_id;
    //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "procedures_name, procedures_type, procedures_dl_name, %s, %s, %s", procedures_name, procedures_type, procedures_dl_name);
    int tmp_control_vars[control_vars.size()+1];
    //tmp_control_vars[0]=this->num_ens_members;
    tmp_control_vars[0]=this->proc_member_id;
    for(int i=1;i<control_vars.size()+1;++i) tmp_control_vars[i]=control_vars[i-1];
    int tmp_decomp_ids[decomp_ids.size()];
    for(int i=0;i<decomp_ids.size();++i) tmp_decomp_ids[i]=decomp_ids[i];
    int tmp_grid_ids[grid_ids.size()];
    for(int i=0;i<grid_ids.size();++i) tmp_grid_ids[i]=grid_ids[i];        
    int tmp_timer_ids[timer_ids.size()+1];
    tmp_timer_ids[0]=this->periodic_timer_id;
    for(int i=1;i<timer_ids.size()+1;++i) tmp_timer_ids[i]=timer_ids[i-1];

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());
    this->external_instance_id = external_procedures_mgr->initialize_external_procedures_inst(instance_name, procedures_name, procedures_type, member_comp_id, procedures_dl_name, 1, control_vars.size()+1, decomp_ids.size(), grid_ids.size(), API_specified_field_insts.size(), timer_ids.size()+1, tmp_control_vars, tmp_decomp_ids, tmp_grid_ids, fields_id, tmp_timer_ids, "Initialize external data assimilation procedures");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());

    External_procedures_inst *external_procedures_inst = external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMM, " Obtain registered external procedures insts");
    external_procedures_inst->get_procedures_import_field_insts(mirror_procedures_import_field_insts);
    external_procedures_inst->get_procedures_export_field_insts(mirror_procedures_export_field_insts);
    //for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++) EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "mirror_procedures_export_field_insts, %s", mirror_procedures_export_field_insts[i]->get_field_name());
    for (int i = 0; i < mirror_procedures_import_field_insts.size(); i ++) {
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_import_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                copy_in_index.push_back(j);
    }
    for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++) {
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_export_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                copy_out_index.push_back(j);
    }
    //mirror_procedures_import_field_insts.clear();
    //mirror_procedures_export_field_insts.clear();	
}


void Ensemble_procedures_inst::do_copy_in()
{   
     for (int i = 0; i < this->mirror_procedures_import_field_insts.size(); i ++) {
        for (int j = 0; j < mirror_API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_import_field_insts[i]->get_field_name(), mirror_API_specified_field_insts[j]->get_field_name()))
                memory_manager->copy_field_data_values(this->mirror_API_specified_field_insts[j], this->API_specified_field_insts[j]);
    }
}


void Ensemble_procedures_inst::do_copy_out()
{   
    for (int i = 0; i < this->mirror_procedures_export_field_insts.size(); i ++) {
        for (int j = 0; j < mirror_API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_export_field_insts[i]->get_field_name(), mirror_API_specified_field_insts[j]->get_field_name()))
                memory_manager->copy_field_data_values(this->API_specified_field_insts[j], this->mirror_API_specified_field_insts[j]);
    }
}


void Ensemble_procedures_inst::do_ensemble_op_initialize()
{
    char tmp_annotation[NAME_STR_SIZE], interface_name[NAME_STR_SIZE];
    int mem_id = this->proc_member_id;
    bool do_ensemble_op_first = true;
    int registered_first_do_ensemble_op_procedure_id;
    for (int i = 1; i < ensemble_procedures_mgr->get_num_registered_ensemble_procedures_insts(); i ++){
        if (ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(i)->get_field_instances_op()->if_do_ensemble_op()){
            do_ensemble_op_first = false;
            registered_first_do_ensemble_op_procedure_id = i;
        }
    }
    if (do_ensemble_op_first) {    
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "do_ensemble_op_first");    
    //注册集合set网格和剖分
    //double value1[this->num_ens_members], temp_value1, *temp_value2;
    //temp_value2 = new double [this->num_ens_members];
    //for (int i = 0; i < this->num_ens_members; i ++) value1[i]=i;
    //transform_datatype_of_arrays((double*)value1, &temp_value1, 1);
    //transform_datatype_of_arrays((double*)value1, temp_value2, this->num_ens_members);
    //int ens_member_virtual_grid_id = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_V1D_Z_GRID_VIA_MODEL, this->set_comp_id, "ens_member_virtual_grid", 1, "num", this->//num_ens_members, temp_value1, temp_value2, NULL, "register the virtual grid with the number of ensemble members for ensemble set");
    int ens_member_virtual_grid_id = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_NORMAL_1D_GRID_NO_DATA, this->set_comp_id, "ens_member_virtual_grid", 0, NULL, this->num_ens_members, 0.0, NULL, NULL, "register the virtual grid with the number of ensemble members for ensemble set");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "this->num_ens_members : %d", this->num_ens_members);
    std::vector<Original_grid_info*> tmp_ensemble_member_grids;
    char set_md_grid_name[NAME_STR_SIZE];
    for (int i = 0; i < this->grid_ids.size(); i ++){
        tmp_ensemble_member_grids.clear();
        tmp_ensemble_member_grids.push_back(original_grid_mgr->promote_ensemble_member_grid_to_set(this->set_comp_id, original_grid_mgr->search_grid_info(grid_ids[i])));
        ensemble_member_grid_ids.push_back(tmp_ensemble_member_grids.back()->get_grid_id());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "register md grid for \"%s\"", tmp_ensemble_member_grids.back()->get_grid_name());
        sprintf(set_md_grid_name, "set_md_grid_%s", tmp_ensemble_member_grids.back()->get_grid_name());
        ensemble_set_grid_ids.push_back(original_grid_mgr->register_md_grid_via_multi_grids(this->set_comp_id, set_md_grid_name, tmp_ensemble_member_grids.back()->get_grid_id(), ens_member_virtual_grid_id, -1, -1, NULL, "register grid with the number of ensemble members for ensemble set"));
        if (tmp_ensemble_member_grids.back()->is_H2D_grid()) {
            this->ensemble_set_decomp_id = decomps_info_mgr->generate_default_parallel_decomp(tmp_ensemble_member_grids.back())->get_decomp_id();
            this->ensemble_set_decomp_ids.push_back(ensemble_set_decomp_id);
        }
        else if(tmp_ensemble_member_grids.back()->is_V3D_grid()){
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "%s", tmp_ensemble_member_grids.back()->get_H2D_sub_CoR_grid()->get_grid_name());
            for(int j = 0; j < this->ensemble_set_decomp_ids.size(); j ++){
                if (decomp_ids[j] != -1){
                EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "%s", original_grid_mgr->search_grid_info(decomps_info_mgr->get_decomp_info(ensemble_set_decomp_ids[j])->get_grid_id())->get_H2D_sub_CoR_grid()->get_grid_name());
                if (words_are_the_same(original_grid_mgr->search_grid_info(decomps_info_mgr->get_decomp_info(ensemble_set_decomp_ids[j])->get_grid_id())->get_H2D_sub_CoR_grid()->get_grid_name(), tmp_ensemble_member_grids.back()->get_H2D_sub_CoR_grid()->get_grid_name())){
                    this->ensemble_set_decomp_ids.push_back(ensemble_set_decomp_ids[j]);
                   break; 
                }
                //else EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens in generating ensemble_set_decomp_ids: the H2D grid must be in front of H3D grid which contains the H2D grid");
                }
            }
        }
        else this->ensemble_set_decomp_ids.push_back(this->decomp_ids[i]);
        //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_member_grid_id: %d", ensemble_member_grid_ids.back());
        //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_set_grid_id: %d", ensemble_set_grid_ids.back());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "j= %d, ensemble_set_decomp_id: %d", i, this->ensemble_set_decomp_ids.back());
    }
    //分别根据集合member和集合set对应的网格和剖分声明field instance（集合member变量共享集合set对应的内存地址）
    int tmp_decomp_id, tmp_set_comp_or_grid_id, tmp_member_comp_or_grid_id;
    char member_grid_name[NAME_STR_SIZE];
    std::vector<Field_mem_info*> tmp_field_insts;
    ensemble_member_field_insts_id = new int *[this->num_ens_members];
    for(int i = 0; i < this->num_ens_members; i++) ensemble_member_field_insts_id[i] = new int [API_specified_field_insts.size()];
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register ensemble set field instance");
    for (int i = 0; i < this->API_specified_field_insts.size(); i ++){
        
        if (API_specified_field_insts[i]->get_grid_id() == -1 || API_specified_field_insts[i]->get_grid_id() == CCPL_NULL_INT ){
            tmp_set_comp_or_grid_id = ens_member_virtual_grid_id;
            tmp_member_comp_or_grid_id = this->set_comp_id;
        }
        for (int j = 0; j < this->ensemble_set_grid_ids.size(); j ++){           
            sprintf(member_grid_name, "set_md_grid_member_grid_%s", API_specified_field_insts[i]->get_grid_name());
            if (words_are_the_same(original_grid_mgr->search_grid_info(ensemble_set_grid_ids[j])->get_grid_name(), member_grid_name)){
                tmp_set_comp_or_grid_id = ensemble_set_grid_ids[j];
                tmp_member_comp_or_grid_id = ensemble_member_grid_ids[j];
                tmp_decomp_id = this->ensemble_set_decomp_ids[j];
            }
            if (words_are_the_same(API_specified_field_insts[i]->get_field_name(), "ccpl_t_global")){
                    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "member_name:%s", member_grid_name);
                    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "%s", original_grid_mgr->search_grid_info(ensemble_set_grid_ids[j])->get_grid_name());
                    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "j= %d, ensemble_set_decomp_id: %d", j, this->ensemble_set_decomp_ids[j]);
                    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "tmp_decomp_id: %d", tmp_decomp_id);
                    }
        }       
        if (API_specified_field_insts[i]->get_decomp_id() == -1 || API_specified_field_insts[i]->get_decomp_id() == CCPL_NULL_INT ) tmp_decomp_id = -1;
        sprintf(tmp_annotation, "Ensemble set: register field instance of %s", API_specified_field_insts[i]->get_field_name());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register ensemble set field instance: %s", API_specified_field_insts[i]->get_field_name());
        ensemble_set_field_insts_before_ens_op.push_back(memory_manager->alloc_mem(API_specified_field_insts[i]->get_field_name(), tmp_decomp_id, tmp_set_comp_or_grid_id, -(this->instance_id), API_specified_field_insts[i]->get_data_type(), API_specified_field_insts[i]->get_unit(), tmp_annotation, false));
        sprintf(tmp_annotation, "Ensemble set: register field instance of %s for ensemble member", API_specified_field_insts[i]->get_field_name());
        ensemble_set_field_insts_one_member.push_back(memory_manager->alloc_mem(API_specified_field_insts[i]->get_field_name(), tmp_decomp_id, tmp_member_comp_or_grid_id, -(this->instance_id+1), API_specified_field_insts[i]->get_data_type(), API_specified_field_insts[i]->get_unit(), tmp_annotation, false));
        ensemble_set_field_insts_one_member.back()->define_field_values(false);
        for (int j = 0; j < this->num_ens_members; j ++){
            tmp_field_insts.clear();
            sprintf(tmp_annotation, "Ensemble set: register field instance of %s for ensemble member_%d", API_specified_field_insts[i]->get_field_name(), j+1);
            tmp_field_insts.push_back(memory_manager->alloc_mem(API_specified_field_insts[i]->get_field_name(), tmp_decomp_id, tmp_member_comp_or_grid_id, -(j+1), API_specified_field_insts[i]->get_data_type(), API_specified_field_insts[i]->get_unit(), tmp_annotation, false));
            tmp_field_insts.back()->reset_mem_buf((char *)(ensemble_set_field_insts_before_ens_op[i]->get_data_buf())+j*tmp_field_insts.back()->get_size_of_field()*get_data_type_size(tmp_field_insts.back()->get_data_type()), true, -1);
            ensemble_member_field_insts_id[j][i] = tmp_field_insts.back()->get_field_instance_id();
        }
    }
    //生成model临时输出接口和各个集合member的临时输入接口
    int tmp_field_update_status_size = API_specified_field_insts.size()+API_specified_field_insts.size()+1;
    int tmp_field_update_status[tmp_field_update_status_size];
    int tmp_ensemble_member_import_interface_id[this->num_ens_members];
    int tmp_model_src_field_ids[API_specified_field_insts.size()];
    int tmp_ensemble_member_dst_field_ids[API_specified_field_insts.size()];
    for (int i = 0; i < API_specified_field_insts.size(); i ++) tmp_model_src_field_ids[i] = API_specified_field_insts[i]->get_field_instance_id();
    sprintf(interface_name, "tmp_model_export_interface_for_member_%d", this->proc_member_id);
    sprintf(tmp_annotation, "register model export interface for ensemble member_%d", this->proc_member_id);
    int tmp_model_export_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 1, API_specified_field_insts.size(), tmp_model_src_field_ids, API_specified_field_insts.size(), timer_ids.front(), 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    //sprintf(tmp_annotation, "execute model export interface for ensemble member_%d", this->proc_member_id);
    //inout_interface_mgr->get_interface(tmp_model_export_interface_id)->execute(true, API_ID_ENSEMBLE_PROC_INST_INIT, tmp_field_update_status, tmp_field_update_status_size, tmp_annotation);
    this->ensemble_set_timer_id = timer_mgr->define_timer(this->set_comp_id, timer_mgr->get_timer(this->timer_ids.back())->get_frequency_unit(), timer_mgr->get_timer(this->timer_ids.back())->get_frequency_count(), timer_mgr->get_timer(this->timer_ids.back())->get_local_lag_count(), timer_mgr->get_timer(this->timer_ids.back())->get_remote_lag_count(), "define timer of ensemble set");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "if this->ensemble_set_timer_id: %d", this->ensemble_set_timer_id);
    for (int i = 0; i < this->num_ens_members; i ++){        
        for (int j = 0; j < API_specified_field_insts.size(); j ++) tmp_ensemble_member_dst_field_ids[j] = ensemble_member_field_insts_id[i][j];
        sprintf(interface_name, "tmp_ensemble_set_import_interface_for_member_%d", i+1);
        sprintf(tmp_annotation, "Ensemble set: register import interface for ensemble member_%d", i+1);
        tmp_ensemble_member_import_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 0, API_specified_field_insts.size(), tmp_ensemble_member_dst_field_ids, API_specified_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
        //sprintf(tmp_annotation, "Ensemble set: execute import interface for ensemble member_%d", i+1);
        //inout_interface_mgr->get_interface(tmp_ensemble_member_import_interface_id[i])->execute(true, API_ID_ENSEMBLE_PROC_INST_INIT, tmp_field_update_status, tmp_field_update_status_size, tmp_annotation);
    }
    //进行model临时输出接口和各个集合member临时输入接口的耦合连接并执行 
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "is_interface_id_legal: %d", inout_interface_mgr->is_interface_id_legal(tmp_model_export_interface_id));
     
    for (int i = 0; i < this->num_ens_members; i ++){
        if (mem_id-1 == i) EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(tmp_ensemble_member_import_interface_id[i]), inout_interface_mgr->get_interface(tmp_model_export_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
        else EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(tmp_ensemble_member_import_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
    }
    
    int num_dst_fields_model, num_dst_fields_ensemble_member;   
    sprintf(tmp_annotation, "execute model export interface for ensemble member_%d", mem_id);
    inout_interface_mgr->execute_interface(tmp_model_export_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, tmp_field_update_status, tmp_field_update_status_size, &num_dst_fields_model, tmp_annotation);
    for (int i = 0; i < this->num_ens_members; i ++){
        sprintf(tmp_annotation, "Ensemble set: execute import interface for ensemble member_%d", i);
        inout_interface_mgr->execute_interface(tmp_ensemble_member_import_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, tmp_field_update_status, tmp_field_update_status_size, &num_dst_fields_ensemble_member, tmp_annotation);
    }
    }
    else{
        for (int i = 1; i <= ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_num_ensemble_member_grid_ids(); i ++) 
            this->ensemble_member_grid_ids.push_back(ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_member_grid_id(i));
        for (int i = 1; i <= ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_num_ensemble_set_grid_ids(); i ++) 
            this->ensemble_set_grid_ids.push_back(ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_set_grid_id(i));
        for (int i = 1; i <= ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_num_ensemble_set_decomp_ids(); i ++) 
            this->ensemble_set_decomp_ids.push_back(ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_set_decomp_id(i));
        for (int i = 1; i <= ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_num_ensemble_set_field_insts_before_ens_op(); i ++) 
            this->ensemble_set_field_insts_before_ens_op.push_back(ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_set_field_insts_before_ens_op(i));
        for (int i = 1; i <= ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_num_ensemble_set_field_insts_one_member(); i ++) 
            this->ensemble_set_field_insts_one_member.push_back(ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_set_field_insts_one_member(i));
        ensemble_member_field_insts_id = new int *[this->num_ens_members];
        for(int i = 0; i < this->num_ens_members; i++) ensemble_member_field_insts_id[i] = new int [API_specified_field_insts.size()];
        for (int i = 0; i < this->API_specified_field_insts.size(); i ++)
            for (int j = 0; j < this->num_ens_members; j ++)
                ensemble_member_field_insts_id[j][i] = ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_member_field_insts_id(i,j);
        this->ensemble_set_timer_id = ensemble_procedures_mgr->get_registered_ensemble_procedures_inst(registered_first_do_ensemble_op_procedure_id)->get_ensemble_set_timer_id();
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "else this->ensemble_set_timer_id: %d", this->ensemble_set_timer_id);
    }
    for(int i=0;i<ensemble_set_field_insts_before_ens_op.size();++i) 
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_set_field_insts_before_ens_op: \"%s\"", ensemble_set_field_insts_before_ens_op[i]->get_field_name());

    //根据xml配置进行集合操作
    do_ensemble_op(API_specified_field_insts, ensemble_set_field_insts_after_ens_op);

    for(int i=0;i<ensemble_set_field_insts_after_ens_op.size();++i) 
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_set_field_insts_after_ens_op: \"%s\"", ensemble_set_field_insts_after_ens_op[i]->get_field_name());

    //通过集合set调用外部过程执行接口
    int tmp_control_vars[control_vars.size()+1];
    tmp_control_vars[0]=this->num_ens_members;
    for(int i=1;i<control_vars.size()+1;++i) tmp_control_vars[i]=control_vars[i-1];
    int tmp_grid_ids[ensemble_member_grid_ids.size()+ensemble_set_grid_ids.size()];
    int tmp_decomp_ids[ensemble_set_decomp_ids.size()+ensemble_set_decomp_ids.size()];
    if (field_instances_op->if_use_set_grids()){
        for(int i=0;i<ensemble_set_grid_ids.size();++i) {
            tmp_grid_ids[i]=ensemble_set_grid_ids[i];
            tmp_decomp_ids[i]=-1;
        }
        for(int i=ensemble_set_grid_ids.size();i<ensemble_member_grid_ids.size()+ensemble_set_grid_ids.size();++i) {
            tmp_grid_ids[i]=ensemble_member_grid_ids[i-ensemble_set_grid_ids.size()];
            tmp_decomp_ids[i]=ensemble_set_decomp_ids[i-ensemble_set_grid_ids.size()];
        }
    }
    else {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "use member grids");
        for(int i=0;i<ensemble_member_grid_ids.size();++i) {
            tmp_grid_ids[i]=ensemble_member_grid_ids[i];
            tmp_decomp_ids[i]=ensemble_set_decomp_ids[i];
        }
        for(int i=ensemble_set_grid_ids.size();i<ensemble_member_grid_ids.size()+ensemble_set_grid_ids.size();++i) {
            tmp_grid_ids[i]=ensemble_set_grid_ids[i-ensemble_set_grid_ids.size()];      
            tmp_decomp_ids[i]=-1;
        }
    }
    int tmp_timer_ids[2];
    tmp_timer_ids[0]=this->periodic_timer_id;
    tmp_timer_ids[1]=this->ensemble_set_timer_id;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "tmp_timer_ids: %d, %d", tmp_timer_ids[0], tmp_timer_ids[1]);
    int tmp_field_ids[ensemble_set_field_insts_after_ens_op.size()];
    for(int i=0;i<ensemble_set_field_insts_after_ens_op.size();++i) {
        tmp_field_ids[i]=ensemble_set_field_insts_after_ens_op[i]->get_field_instance_id();
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_set_field_insts_after_ens_op: \"%s\"", ensemble_set_field_insts_after_ens_op[i]->get_field_name());
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());
    this->external_instance_id = external_procedures_mgr->initialize_external_procedures_inst(instance_name, procedures_name, procedures_type, set_comp_id, procedures_dl_name, 1, control_vars.size()+1, ensemble_set_decomp_ids.size()+ensemble_set_decomp_ids.size(), ensemble_member_grid_ids.size()+ensemble_set_grid_ids.size(), ensemble_set_field_insts_after_ens_op.size(), 2, tmp_control_vars, tmp_decomp_ids, tmp_grid_ids, tmp_field_ids, tmp_timer_ids, "Initialize external data assimilation procedures");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());

    //获取外部过程field instance输入输出参数列表并生成model和集合member的输入输出列表
    External_procedures_inst *external_procedures_inst = external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMM, " Obtain registered external procedures insts");
    external_procedures_inst->get_procedures_import_field_insts(mirror_procedures_import_field_insts);
    external_procedures_inst->get_procedures_export_field_insts(mirror_procedures_export_field_insts);
    //for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++) EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "mirror_procedures_export_field_insts, %s", mirror_procedures_export_field_insts[i]->get_field_name());
    model_export_field_insts_id = new int [mirror_procedures_import_field_insts.size()];
    ensemble_member_import_field_insts_id = new int *[this->num_ens_members];
    for(int i = 0; i < this->num_ens_members; i++) ensemble_member_import_field_insts_id[i] = new int [mirror_procedures_import_field_insts.size()];
    for (int i = 0; i < mirror_procedures_import_field_insts.size(); i ++) {
        for (int j = 0; j < API_specified_field_insts.size(); j ++){
            if (words_are_the_same(mirror_procedures_import_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name())){
                model_export_field_insts_id[i] = API_specified_field_insts[j]->get_field_instance_id();
                for (int k = 0; k < this->num_ens_members; k++) ensemble_member_import_field_insts_id[k][i] = ensemble_member_field_insts_id[k][j];
            }
        }
    }
    model_import_field_insts_id = new int [mirror_procedures_export_field_insts.size()];
    ensemble_member_export_field_insts_id = new int *[this->num_ens_members];
    for(int i = 0; i < this->num_ens_members; i++) ensemble_member_export_field_insts_id[i] = new int [mirror_procedures_export_field_insts.size()];
    for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++) {
        for (int j = 0; j < API_specified_field_insts.size(); j ++){
            if (words_are_the_same(mirror_procedures_export_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name())){
                model_import_field_insts_id[i] = API_specified_field_insts[j]->get_field_instance_id();
                for (int k = 0; k < this->num_ens_members; k++) ensemble_member_export_field_insts_id[k][i] = ensemble_member_field_insts_id[k][j];
            }
        }
    }
    //生成新的model输入输出接口和各个集合member的输入输出接口
    sprintf(interface_name, "%s.model_export_interface_for_member_%d", this->instance_name, this->proc_member_id);
    sprintf(tmp_annotation, "register model export interface for ensemble member_%d", this->proc_member_id);
    this->model_export_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 1, mirror_procedures_import_field_insts.size(), model_export_field_insts_id, mirror_procedures_import_field_insts.size(), timer_ids.front(), 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    if (!mirror_procedures_export_field_insts.empty()){
        sprintf(interface_name, "%s.model_import_interface_for_member_%d", this->instance_name, this->proc_member_id);
        sprintf(tmp_annotation, "register model import interface for ensemble member_%d", this->proc_member_id);
        this->model_import_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 0, mirror_procedures_export_field_insts.size(), model_import_field_insts_id, mirror_procedures_export_field_insts.size(), timer_ids.front(), 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    }

    int tmp_ensemble_member_import_field_ids[mirror_procedures_import_field_insts.size()];
    ensemble_member_import_interface_id = new int [this->num_ens_members];
    int tmp_ensemble_member_export_field_ids[mirror_procedures_export_field_insts.size()];
    ensemble_member_export_interface_id = new int [this->num_ens_members];
    for (int i = 0; i < this->num_ens_members; i ++){        
        for (int j = 0; j < mirror_procedures_import_field_insts.size(); j ++) tmp_ensemble_member_import_field_ids[j] = ensemble_member_import_field_insts_id[i][j];
        sprintf(interface_name, "%s.ensemble_set_import_interface_for_member_%d", this->instance_name, i+1);
        sprintf(tmp_annotation, "Ensemble set: register import interface for ensemble member_%d", i+1);
        ensemble_member_import_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 0, mirror_procedures_import_field_insts.size(), tmp_ensemble_member_import_field_ids, mirror_procedures_import_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
        if (!mirror_procedures_export_field_insts.empty()){
            for (int j = 0; j < mirror_procedures_export_field_insts.size(); j ++) tmp_ensemble_member_export_field_ids[j] = ensemble_member_export_field_insts_id[i][j];
            sprintf(interface_name, "%s.ensemble_set_export_interface_for_member_%d", this->instance_name, i+1);
            sprintf(tmp_annotation, "Ensemble set: register export interface for ensemble member_%d", i+1);
            ensemble_member_export_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 1, mirror_procedures_export_field_insts.size(), tmp_ensemble_member_export_field_ids, mirror_procedures_export_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER); 
        }  
    }
    //进行model输入输出接口和各个集合member输入输出接口的耦合连接
    for (int i = 0; i < this->num_ens_members; i ++){
        if (mem_id-1 == i) {
            EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_import_interface_id[i]), inout_interface_mgr->get_interface(this->model_export_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
            if (!mirror_procedures_export_field_insts.empty()) 
                EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i]), inout_interface_mgr->get_interface(this->model_import_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
        }
        else {
            EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_import_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
            if (!mirror_procedures_export_field_insts.empty()) 
                EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_op_initialize");
        }
    }
}


void Ensemble_procedures_inst::execute_config_script(int comp_id, const char *file_name, const char *str_para0, const char *str_para1, const char *str_para2, const char *str_para3, const char *str_para4, const char *str_para5, const char *annotation)
{
    int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "Ensemble_procedures_inst::execute_config_script");
    //synchronize_comp_processes_for_API(comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"), "synchorization before executing config script", "Ensemble_procedures_inst::execute_config_script");
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"));
    if ( local_proc_id == 0 ) {

        char full_file_name[NAME_STR_SIZE*16], working_directory[NAME_STR_SIZE*16], full_command[NAME_STR_SIZE*32], tmp_full_command[NAME_STR_SIZE*32];
        const char *str_paras[6];
    
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "execute_config_script for \"%s\" begins at the model code with the annotation \"%s\"", file_name, annotation);
        realpath(file_name, full_file_name);
        getcwd(working_directory, NAME_STR_SIZE*16);
        EXECUTION_REPORT(REPORT_ERROR, -1, does_file_exist(file_name), "ERROR happens when executing the execute_config_script at the model code with the annotation \"%s\" for the script \"%s\": this script (the absolute file name is \"%s\") does not exist. Please verify (please note that the current working directory is \"%s\").", annotation, file_name, full_file_name,    working_directory);
        str_paras[0] = str_para0;
        str_paras[1] = str_para1;
        str_paras[2] = str_para2;
        str_paras[3] = str_para3;
        str_paras[4] = str_para4;
        str_paras[5] = str_para5;
        sprintf(full_command, "\"%s\"", full_file_name);
        for (int i = 0; i < 5; i ++)
            if (strlen(str_paras[i]) > 0) {
                strcpy(tmp_full_command, full_command);
                sprintf(full_command, "%s \"%s\"", tmp_full_command, str_paras[i]);
            }
        EXECUTION_REPORT(REPORT_ERROR, -1, system(full_command) != -1, "ERROR happens when executing the execute_config_script at the model code with the annotation \"%s\" for the script \"%s\"    (\"%s\"): fail to execute this script. Please verify.", annotation, file_name, full_file_name);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "execute_config_script for \"%s\" ends", file_name);
    }
    //synchronize_comp_processes_for_API(comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"), "synchorization after executing config script", "Ensemble_procedures_inst::execute_config_script");
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"));
}


void Ensemble_procedures_inst::run(bool do_advance_time, int chunk_index, const char *annotation)
{
    check_for_ccpl_managers_allocated(API_ID_TIME_MGT_IS_TIMER_ON, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to Ensemble_procedures_inst::run: %s@%s", this->instance_name, this->get_procedures_name());
    
    if (timer_mgr->is_timer_on(this->periodic_timer_id, annotation)){
        char full_date[15], date[9], hours[3], minutes[3], seconds[3];
        sprintf(date, "%08d", components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_date());
        std::strncpy(date, date, 8);
        sprintf(hours, "%02d", (components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second())/3600);
        std::strncpy(hours, hours, 2);
        sprintf(minutes, "%02d", ((components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second())%3600)/60);
        std::strncpy(minutes, minutes, 2);
        sprintf(seconds, "%02d", ((components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second())%3600)%60);
        std::strncpy(seconds, seconds, 2);
        sprintf(full_date, "%08s%02s%02s%02s", date, hours, minutes, seconds);
        std::strncpy(full_date, full_date, 14);
       
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "date: %08s", date);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "length_date: %d", strlen(date));
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "hours: %02s", hours);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "minutes: %02s", minutes);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "seconds: %02s", seconds);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "current full date: %14s", full_date);
        if (this->get_field_instances_op()->if_do_field_instances_operation() && !this->get_field_instances_op()->if_do_ensemble_op()){
            this->do_copy_in();  
            if (this->before_instance_script != NULL ) execute_config_script(this->set_comp_id, this->before_instance_script, full_date, "", "", "", "", "", "before instance configuration script");
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            //synchronize_comp_processes_for_API(this->member_comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"), "synchorization before running external procedures", "Ensemble_procedures_inst::run");
            check_for_ccpl_managers_allocated(API_ID_EXTERNAL_PROC_INST_RUN, annotation);
            external_procedures_mgr->get_procedures_inst(this->get_external_instance_id(), API_ID_EXTERNAL_PROC_INST_RUN, annotation)->run(chunk_index, annotation);
            //synchronize_comp_processes_for_API(this->member_comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"), "synchorization after running external procedures", "Ensemble_procedures_inst::run");
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            if (this->after_instance_script != NULL ) execute_config_script(this->set_comp_id, this->after_instance_script, full_date, "", "", "", "", "", "after instance configuration script");
            this->do_copy_out();
        }
        if (this->get_field_instances_op()->if_do_field_instances_operation() && this->get_field_instances_op()->if_do_ensemble_op()){
            char tmp_annotation[NAME_STR_SIZE];
            int num_dst_fields_model, num_dst_fields_ensemble_member;  
            int model_export_field_update_status[mirror_procedures_import_field_insts.size()+mirror_procedures_export_field_insts.size()+1]; 
            int mem_id = this->proc_member_id;
            sprintf(tmp_annotation, "execute model export interface for ensemble member_%d", this->proc_member_id);
            inout_interface_mgr->execute_interface(this->model_export_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size()+mirror_procedures_export_field_insts.size()+1, &num_dst_fields_model, tmp_annotation);
            for (int i = 0; i < this->num_ens_members; i ++){
                sprintf(tmp_annotation, "Ensemble set: execute import interface for ensemble member_%d", i);
                inout_interface_mgr->execute_interface(ensemble_member_import_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size()+mirror_procedures_export_field_insts.size()+1, &num_dst_fields_ensemble_member, tmp_annotation);
            }
            ensemble_set_field_insts_after_ens_op.clear();
            do_ensemble_op(mirror_procedures_import_field_insts, ensemble_set_field_insts_after_ens_op);
            if (this->before_instance_script != NULL ) execute_config_script(this->set_comp_id, this->before_instance_script, full_date, "", "", "", "", "", "before instance configuration script");
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            //synchronize_comp_processes_for_API(this->set_comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->set_comp_id, "Ensemble_procedures_inst::run"), "synchorization before running external procedures", "Ensemble_procedures_inst::run");
            check_for_ccpl_managers_allocated(API_ID_EXTERNAL_PROC_INST_RUN, annotation);
            external_procedures_mgr->get_procedures_inst(this->get_external_instance_id(), API_ID_EXTERNAL_PROC_INST_RUN, annotation)->run(chunk_index, annotation);
            //synchronize_comp_processes_for_API(this->set_comp_id, API_ID_COMP_MGT_END_COMP_REG, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->set_comp_id, "Ensemble_procedures_inst::run"), "synchorization after running external procedures", "Ensemble_procedures_inst::run");
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            if (this->after_instance_script != NULL ) execute_config_script(this->set_comp_id, this->after_instance_script, full_date, "", "", "", "", "", "after instance configuration script");
            if (!mirror_procedures_export_field_insts.empty()){
                int model_import_field_update_status[mirror_procedures_import_field_insts.size()+mirror_procedures_export_field_insts.size()+1];
                for (int i = 0; i < this->num_ens_members; i ++){
                    sprintf(tmp_annotation, "Ensemble set: execute export interface for ensemble member_%d", i);
                    inout_interface_mgr->execute_interface(ensemble_member_export_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_import_field_update_status, mirror_procedures_export_field_insts.size()+mirror_procedures_export_field_insts.size()+1, &num_dst_fields_ensemble_member, tmp_annotation);
                    if (mem_id-1 == i){
                        sprintf(tmp_annotation, "execute model import interface for ensemble member_%d", this->proc_member_id);
                        inout_interface_mgr->execute_interface(this->model_import_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_import_field_update_status, mirror_procedures_import_field_insts.size()+mirror_procedures_export_field_insts.size()+1, &num_dst_fields_model, tmp_annotation);
                    }
                }
            }
        }
    }
    if (do_advance_time){
        components_time_mgrs->advance_component_time(this->member_comp_id, "member advance time");
        EXECUTION_REPORT(REPORT_PROGRESS, this->member_comp_id, true, "Component model \"%s\" advance time at the model code with the annotation \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(this->member_comp_id,true,"")->get_full_name(), "member advance time");
        components_time_mgrs->advance_component_time(this->set_comp_id, "ensemble set advance time");
        EXECUTION_REPORT(REPORT_PROGRESS, this->set_comp_id, true, "Component model \"%s\" advance time at the model code with the annotation \"%s\"", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(this->set_comp_id,true,"")->get_full_name(), "ensemble set advance time");
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "date: %d",  components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_date());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "seconds: %d",  components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "date: %d",  components_time_mgrs->get_time_mgr(this->set_comp_id)->get_current_date());
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "seconds: %d",  components_time_mgrs->get_time_mgr(this->set_comp_id)->get_current_second());
    }
    
}

Ensemble_procedures_mgt::Ensemble_procedures_mgt()
{
}


Ensemble_procedures_mgt::~Ensemble_procedures_mgt()
{
	for (int i = 0; i < registered_ensemble_procedures_insts.size(); i ++)
		delete registered_ensemble_procedures_insts[i];
}

int Ensemble_procedures_mgt::initialize_ensemble_procedures_inst(const char *inst_name, int set_comp_id, int member_comp_id, int size_field_inst, int size_grids, int size_decomps, int size_timers, int size_controls,  
                                                                 const int *field_inst_ids, const int *grid_ids, const int *decomp_ids, const int *timer_ids, const int *control_vars, const char *annotation)
{
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start Ensemble_procedures_mgt::initialize_ensemble_procedures_inst");
    this->registered_ensemble_procedures_insts_run_count = 0;
	for (int i = 0; i < registered_ensemble_procedures_insts.size(); i ++)
		if (registered_ensemble_procedures_insts[i]->get_set_comp_id() == set_comp_id && words_are_the_same(registered_ensemble_procedures_insts[i]->get_instance_name(), inst_name)) 
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, set_comp_id, false, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\" at the model code with the annotation \"%s\": annother instance with the same name has been initialized before (at the model code with the annotation \"%s\"). Please verify.", inst_name, annotation, annotation_mgr->get_annotation(registered_ensemble_procedures_insts[i]->get_instance_id(), "registering an instance of ensemble procedures")); 
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "registered_ensemble_procedures_insts.size is \"%d\"",registered_ensemble_procedures_insts.size());
	int instance_id = TYPE_ENS_PROCEDURE_PREFIX | (registered_ensemble_procedures_insts.size());
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to Ensemble_procedures_inst");
	Ensemble_procedures_inst *procedures_inst = new Ensemble_procedures_inst(instance_id, inst_name, set_comp_id, member_comp_id, size_field_inst, size_grids, size_decomps, size_timers, size_controls, 
																			 field_inst_ids, grid_ids, decomp_ids, timer_ids, control_vars, annotation);	
	return instance_id;
}

Ensemble_procedures_inst *Ensemble_procedures_mgt::get_procedures_inst(int instance_id, int API_id, const char *annotation)
{
	char API_label[NAME_STR_SIZE]; 
	int instance_index;

	
	get_API_hint(-1, API_id, API_label);
	instance_index = GET_ENS_PROCEDURES_INST_INDEX(instance_id);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (instance_id&TYPE_ID_PREFIX_MASK) == TYPE_ENS_PROCEDURE_PREFIX && instance_index < registered_ensemble_procedures_insts.size(), "ERROR happens when calling the API \"%s\": the given \"instance_id\" (%x) is not a legal ID of an instance of ensemble procedures. Please verify the model code with the anotation \"%s\"", API_label, instance_id, annotation);

	return registered_ensemble_procedures_insts[instance_index];
}


MPI_Comm Ensemble_procedures_mgt::get_instance_local_comm(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMM, annotation)->get_local_comm();
}


int Ensemble_procedures_mgt::get_instance_member_comp_id(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMP, annotation)->get_member_comp_id();
}

int Ensemble_procedures_mgt::get_instance_set_comp_id(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMP, annotation)->get_set_comp_id();
}

int Ensemble_procedures_mgt::get_instance_num_grid_decomps(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_NUM_GRID_DECOMPS, annotation)->get_num_grid_decomps();
}


int Ensemble_procedures_mgt::get_instance_grid_id(int instance_id, int grid_index, const char *annotation)
{
	Ensemble_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_GRID_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), grid_index >= 1, "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of grid_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), grid_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), grid_index <= procedures_inst->get_num_grid_decomps(), "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_grid_ID\" regarding the ensemble procedures instance \"%s\": the value of grid_index (%d) cannot be larger than %d (the number of given grids when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), grid_index, procedures_inst->get_num_grid_decomps(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_grid_id(grid_index);
}


int Ensemble_procedures_mgt::get_instance_decomp_id(int instance_id, int decomp_index, const char *annotation)
{
	Ensemble_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_DECOMP_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), decomp_index >= 1, "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of decomp_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), decomp_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), decomp_index <= procedures_inst->get_num_grid_decomps(), "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_decomp_ID\" regarding the ensemble procedures instance \"%s\": the value of decomp_index (%d) cannot be larger than %d (the number of given parallel decompositions when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), decomp_index, procedures_inst->get_num_grid_decomps(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_decomp_id(decomp_index);
}


int Ensemble_procedures_mgt::get_instance_num_control_vars(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_NUM_CONTROLS, annotation)->get_num_control_vars();
}


int Ensemble_procedures_mgt::get_instance_control_var(int instance_id, int control_var_index, const char *annotation)
{
	Ensemble_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_CONTROL_VAR, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), control_var_index >= 1, "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of control_var_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), control_var_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), control_var_index <= procedures_inst->get_num_control_vars(), "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_control_var\" regarding the ensemble procedures instance \"%s\": the value of control_var_index (%d) cannot be larger than %d (the number of given control variables when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), control_var_index, procedures_inst->get_num_control_vars(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_control_var(control_var_index);
}


int Ensemble_procedures_mgt::get_instance_num_timers(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_NUM_TIMERS, annotation)->get_num_timers();
}


int Ensemble_procedures_mgt::get_instance_timer_id(int instance_id, int timer_index, const char *annotation)
{
	Ensemble_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_TIMER_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), timer_index >= 1, "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of timer_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), timer_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), timer_index <= procedures_inst->get_num_timers(), "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_timer_ID\" regarding the ensemble procedures instance \"%s\": the value of timer_index (%d) cannot be larger than %d (the number of given timers when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), timer_index, procedures_inst->get_num_timers(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_timer_id(timer_index);
}


int Ensemble_procedures_mgt::get_instance_num_fields(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_NUM_FIELDS, annotation)->get_num_specified_field_instances();
}


int Ensemble_procedures_mgt::get_instance_field_id(int instance_id, int field_index, const char *annotation)
{
	Ensemble_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_ENSEMBLE_PROC_INST_GET_FIELD_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), field_index >= 1, "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of field_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), field_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_member_comp_id(), field_index <= procedures_inst->get_num_specified_field_instances(), "ERROR happens when calling the API \"CCPL_ensemble_procedures_para_get_field_ID\" regarding the ensemble procedures instance \"%s\": the value of field_index (%d) cannot be larger than %d (the number of given fields when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), field_index, procedures_inst->get_num_specified_field_instances(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_specified_field_instance_id(field_index);
}

void Ensemble_procedures_mgt::set_registered_ensemble_procedures_insts_run_count()
{
    this->registered_ensemble_procedures_insts_run_count = this->registered_ensemble_procedures_insts_run_count+1;
}
