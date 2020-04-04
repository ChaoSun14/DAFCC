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
#include "external_procedures_mgt.h"



External_procedure_fields_names_mappers::External_procedure_fields_names_mappers(int comp_id, TiXmlElement *mappers_element, const char *XML_file_name)
{
	std::pair<const char *, std::vector<std::pair<const char *, const char*> > > fields_names_mapper;
	

	if (!is_XML_setting_on(comp_id, mappers_element, XML_file_name, "the status of \"fields_name_mappers\"", "external procedures configuration"))
		return;
		
	this->XML_file_name = strdup(XML_file_name);
	for (TiXmlNode *mapper_node = mappers_element->FirstChild(); mapper_node != NULL; mapper_node = mapper_node->NextSibling()) {
		TiXmlElement *mapper_element = mapper_node->ToElement();
		EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(mapper_element->Value(),"fields_name_mapper"), "Detect a wrong name (\"%s\") of an XML element in the XML configuration file \"%s\": the right name should be \"fields_name_mapper\". Please verify the XML file arround the line number %d.", mapper_element->Value(), XML_file_name, mapper_element->Row());
		if (!is_XML_setting_on(comp_id, mapper_element, XML_file_name, "the status of \"fields_name_mapper\"", "external procedures configuration"))
			continue;
		fields_names_mapper.first = NULL;
		fields_names_mapper.second.clear();
		const char *mapper_name = get_XML_attribute(comp_id, -1, mapper_element, "mapper_name", XML_file_name, line_number, "the name of the fields name mapper", "external procedures configuration", true);
		fields_names_mapper.first = strdup(mapper_name);		
		for (TiXmlNode *entry_node = mapper_element->FirstChild(); entry_node != NULL; entry_node = entry_node->NextSibling()) {
			TiXmlElement *entry_element = entry_node->ToElement();
			const char *model_field_name = get_XML_attribute(comp_id, -1, entry_element, "model_field_name", XML_file_name, line_number, "the field name in the model", "external procedures configuration", true);
			const char *procedures_field_name = get_XML_attribute(comp_id, -1, entry_element, "procedures_field_name", XML_file_name, line_number, "the field name in the external procedures", "external procedures configuration", true);
			for (int i = 0; i < (fields_names_mapper.second).size(); i ++) 
				if (words_are_the_same((fields_names_mapper.second)[i].second, procedures_field_name))
					if (words_are_the_same((fields_names_mapper.second)[i].first, model_field_name)) 
						continue;
					else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the number %d field name mapper (\"%s\" <==> \"%s\") is the same or has conflicts with the number %d field name mapper (\"%s\" <==> \"%s\"). Please verify.", entry_element->Row(), XML_file_name, i+1, (fields_names_mapper.second)[i].first, (fields_names_mapper.second)[i].second, (fields_names_mapper.second).size()+1, model_field_name, procedures_field_name);
			std::pair<const char *, const char *> field_names_mapper;
			field_names_mapper.first = strdup(model_field_name);
			field_names_mapper.second = strdup(procedures_field_name);
			(fields_names_mapper.second).push_back(field_names_mapper);			
		}
		fields_names_mappers.push_back(fields_names_mapper);
	}
}


External_procedure_fields_names_mappers::~External_procedure_fields_names_mappers()
{
	delete [] XML_file_name;
	
	for (int i = 0; i < fields_names_mappers.size(); i ++) {
		for (int j = 0; j < fields_names_mappers[i].second.size(); j ++) {
			delete [] (fields_names_mappers[i].second)[j].first;
			delete [] (fields_names_mappers[i].second)[j].second;
		}
		delete [] fields_names_mappers[i].first;
	}
}


std::vector<std::pair<const char *, const char*> > *External_procedure_fields_names_mappers::search_fields_name_mapper(const char *mapper_name)
{
	for (int i = 0; i < fields_names_mappers.size(); i ++) {
		if (words_are_the_same(fields_names_mappers[i].first, mapper_name))
			return &(fields_names_mappers[i].second);
	}

	return NULL;
}


External_procedure_inst::External_procedure_inst(External_procedures_inst *external_procedures_inst, int instance_id, const char *inst_name, const char *procedure_name, int host_comp_id, const char *dl_name, std::vector<std::pair<const char *, const char*> > *fields_names_mapper, const char *annotation)
{
	char tmp_string[NAME_STR_SIZE];
	void *current_dl_handler;


	this->instance_id = instance_id;
	this->procedure_name = strdup(procedure_name);
	local_dl_handler = external_procedures_mgr->get_or_open_dl_handler(host_comp_id, inst_name, procedure_name, dl_name, (char**)(&(this->dl_name)), annotation);
	this->external_procedures_inst = external_procedures_inst;
	this->host_comp_id = external_procedures_inst->host_comp_id;
	this->fields_names_mapper = fields_names_mapper;
	this->require_coupling_parameters = false;

	sprintf(tmp_string, "%s_ccpl_init", procedure_name);
	*((void**)(&procedure_init)) = dlsym(local_dl_handler, tmp_string);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, procedure_init != NULL, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the initialization subroutine \"%s\" does not exist in the dynamic-link library (%s). Please verify at the model code with the annotation \"%s\" or the dynamic-link library.", inst_name, procedure_name, tmp_string, dlerror(), annotation);
	sprintf(tmp_string, "%s_ccpl_run", procedure_name);
	*((void**)(&procedure_run)) = dlsym(local_dl_handler, tmp_string);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, procedure_run != NULL, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the running subroutine \"%s\" does not exist in the dynamic-link library (%s). Please verify at the model code with the annotation \"%s\" or the dynamic-link library.", inst_name, procedure_name, tmp_string, dlerror(), annotation);
	sprintf(tmp_string, "%s_ccpl_finalize", procedure_name);
	*((void**)(&procedure_finalize)) = dlsym(local_dl_handler, tmp_string);

	this->external_procedures_inst->add_external_procedure_inst(this);
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Before executing the external procedures of initializing \"%s\" corresponding to the instance id %d", procedure_name, instance_id);
	procedure_init(&instance_id);
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "After executing the external procedures of initializing \"%s\" corresponding to the instance id %d", procedure_name, instance_id);
	inout_interface_before_run = NULL;
	inout_interface_after_run = NULL;
	field_update_status = new int [procedure_import_field_insts.size()+procedure_export_field_insts.size()+1];
	for (int i = 0; i < declared_internal_fields.size(); i ++) 
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, declared_internal_fields[i].first->get_num_chunks() == 0 || declared_internal_fields[i].first->get_num_chunks() == declared_internal_fields[i].second.size(), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\" at the model code with the annotation \"%s\": field \"%s\" that has been declared as a parameter field has %d chunks while only %d chunks has been declared in the code of the external procedure. Please verify. ", inst_name, procedure_name, annotation, declared_internal_fields[i].first->get_field_name(), declared_internal_fields[i].first->get_num_chunks(), declared_internal_fields[i].second.size());
	if (procedure_import_field_insts.size() > 0) {
		int *src_field_ids = new int [procedure_import_field_insts.size()];
		int *dst_field_ids = new int [procedure_import_field_insts.size()];
		char interface_name[NAME_STR_SIZE];
		sprintf(interface_name, "first parameter interface for procedure %s of procedures instance %s for \"%s\" with index %d", procedure_name, inst_name, comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_full_name(), GET_INST_PROCEDURE_INDEX(instance_id));
		for (int i = 0; i < procedure_import_field_insts.size(); i ++) {
			src_field_ids[i] = model_export_field_insts[i]->get_field_instance_id();
			dst_field_ids[i] = procedure_import_field_insts[i]->get_field_instance_id();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(model_export_field_insts[i]->get_field_name(), procedure_import_field_insts[i]->get_field_name()), "software error in External_procedure_inst::External_procedure_inst");
		}
		int interface_id = inout_interface_mgr->register_normal_remap_interface(interface_name, procedure_import_field_insts.size(), src_field_ids, dst_field_ids, timer_mgr->define_pseudo_timer(host_comp_id), USING_INSTANTANEOUS_VALUE, procedure_import_field_insts.size(), procedure_import_field_insts.size(), API_ID_EXTERNAL_PROC_INST_INIT, annotation);
		inout_interface_before_run = inout_interface_mgr->get_interface(interface_id);
		inout_interface_before_run->execute(true, API_ID_EXTERNAL_PROC_INST_INIT, field_update_status, procedure_import_field_insts.size()+procedure_export_field_insts.size()+1, annotation);
		delete [] src_field_ids;
		delete [] dst_field_ids;
	}
	if (procedure_export_field_insts.size() > 0) {
		int *src_field_ids = new int [procedure_export_field_insts.size()];
		int *dst_field_ids = new int [procedure_export_field_insts.size()];
		char interface_name[NAME_STR_SIZE];
		sprintf(interface_name, "second parameter interface for procedure %s of procedures instance %s for \"%s\" with index %d", procedure_name, inst_name, comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_full_name(), GET_INST_PROCEDURE_INDEX(instance_id));
		for (int i = 0; i < procedure_export_field_insts.size(); i ++) {
			src_field_ids[i] = procedure_export_field_insts[i]->get_field_instance_id();
			dst_field_ids[i] = model_import_field_insts[i]->get_field_instance_id();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(model_import_field_insts[i]->get_field_name(), procedure_export_field_insts[i]->get_field_name()), "software error in External_procedure_inst::External_procedure_inst");
		}
		int interface_id = inout_interface_mgr->register_normal_remap_interface(interface_name, procedure_export_field_insts.size(), src_field_ids, dst_field_ids, timer_mgr->define_pseudo_timer(host_comp_id), USING_INSTANTANEOUS_VALUE, procedure_export_field_insts.size(), procedure_export_field_insts.size(), API_ID_EXTERNAL_PROC_INST_INIT, annotation);
		inout_interface_after_run = inout_interface_mgr->get_interface(interface_id);
		delete [] src_field_ids;
		delete [] dst_field_ids;
	}
}


External_procedure_inst::~External_procedure_inst()
{
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Before calling finalization subroutine of an instance of the external procedure \"%s\"", procedure_name);
	if (procedure_finalize != NULL)
		procedure_finalize(&instance_id);
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Finish finalization subroutine of an instance of the external procedure \"%s\"", procedure_name);
	delete [] dl_name;
	delete [] procedure_name;
	delete [] field_update_status;
	if (inout_interface_before_run != NULL)
		inout_interface_mgr->finalize_inout_interface(inout_interface_before_run->get_interface_id());
	if (inout_interface_after_run != NULL)
		inout_interface_mgr->finalize_inout_interface(inout_interface_after_run->get_interface_id());
}


int External_procedure_inst::declare_para_field(void **data_ptr, const char *field_name, const char *data_type, int type_inout, int decomp_id, int grid_id, int given_num_dims, int *dims_size, int data_size, const char *field_unit, const char *annotation)
{
	std::vector<Field_mem_info*> comp_field_insts;
	char tmp_str1[NAME_STR_SIZE*32], tmp_str2[NAME_STR_SIZE*32], tag[NAME_STR_SIZE], decomp_name[NAME_STR_SIZE], grid_name[NAME_STR_SIZE];
	Decomp_info *decomp_info = NULL;
	Original_grid_info *original_grid = NULL;
	int existing_declared_field_index = -1, comp_or_grid_id, specified_size = 1;
	std::pair<Field_mem_info*, std::vector<std::vector<int> > > declared_internal_field;
	std::vector<int> vec_dims_size;
	int num_existing_chunks;
	bool use_the_same_field_inst;
	MPI_Comm local_comm = external_procedures_inst->local_comm;


	for (int i = 0; i < declared_internal_fields.size(); i ++)
		if (words_are_the_same(declared_internal_fields[i].first->get_field_name(), lookup_model_field_name(field_name))) {
			existing_declared_field_index = i;
			break;
		}

	if (existing_declared_field_index == -1) {
	    synchronize_comp_processes_for_API(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", annotation);
    	check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", field_name, "field_name", annotation);
	    check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", data_type, "implicit data type", annotation);
		check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", type_inout, "type_inout", annotation);	
		check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", given_num_dims, "implict number dimensions of the given pointer", annotation);
		check_API_parameter_decomp(host_comp_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", decomp_id, "decomp_id", annotation);
		int comp_or_grid_id = grid_id;
		if (grid_id == -1)
			comp_or_grid_id = host_comp_id;
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid_mgr->is_grid_id_legal(grid_id), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\" (model field \"%s\" after name mapping): the given \"grid_id\" is not a legal ID of grids. Please verify the model code with the annotation of \"%s\".", field_name, lookup_model_field_name(field_name), annotation);
		check_API_parameter_comp_or_grid(comp_or_grid_id, API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD, local_comm, "declaring a parameter field of external procedures", "comp_or_grid_id", annotation);
	}

	for (int i = 0; i < given_num_dims; i ++)
		specified_size *= dims_size[i];
	if (data_size != -1)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, data_size == specified_size, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\" (model field \"%s\" after name mapping): The size (%d) of the parameter \"data_pointer\" is different from the data size (%d) calculated from of the input parameter \"dims_size\". Please verify the model code with the annotation of \"%s\"" , field_name, lookup_model_field_name(field_name), data_size, specified_size, annotation);
	if (grid_id != -1) {
		original_grid = original_grid_mgr->get_original_grid(grid_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\": the input parameter \"grid_id\" specifies a grid in a different component model (model full name is \"%s\", grid name is \"%s\"). Please verify the model code with the annotation of \"%s\".", field_name, original_grid->get_comp_full_name(), original_grid->get_grid_name(), annotation);
	}
	if (decomp_id != -1) {
		decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, decomp_info->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\" (model field \"%s\" after name mapping): the input parameter \"decomp_id\" specifies a parallel decomposition in a different component model (model full name is \"%s\", parallel decomposition name is \"%s\"). Please verify the model code with the annotation of \"%s\".", field_name, lookup_model_field_name(field_name), comp_comm_group_mgt_mgr->search_global_node(decomp_info->get_comp_id())->get_full_name(), decomp_info->get_decomp_name(), annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid != NULL && original_grid_mgr->get_original_CoR_grid(decomp_info->get_grid_id())->is_subset_of_grid(original_grid->get_original_CoR_grid()), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\" (model field \"%s\" after name mapping): the grid (\"%s\") corresponding to the parallel decomposition \"decomp_id\" is not a subgrid of the grid (\"%s\") corresponding \"grid_id\". Please check the model code with the annotation \"%s\"", field_name, lookup_model_field_name(field_name), decomp_info->get_grid_name(), original_grid->get_grid_name(), annotation);
	}
	sprintf(tag, "declear field %s", lookup_model_field_name(field_name));
	if (decomp_info == NULL || decomp_info->get_num_chunks() == 0) {
		if (existing_declared_field_index != -1)
			EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": this field has been declared before at the model code with the annotation \"%s\" (an instance of an external procedure can declare a field at most once). Please verify. ", procedure_name, field_name, lookup_model_field_name(field_name), annotation, annotation_mgr->get_annotation(instance_id, tag));
		annotation_mgr->add_annotation(instance_id, tag, annotation);
	}
	else {
		if (existing_declared_field_index == -1)
			annotation_mgr->add_annotation(instance_id, tag, annotation);
		else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, declared_internal_fields[existing_declared_field_index].second.size() < decomp_info->get_num_chunks(), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": cannot declare one more chunk of this field because all of the %d chunks (the corresponding parallel decomposition \"%s\" has only %d chunks) have already been declared in the code of the external procedure. Please verify. ", procedure_name, field_name, lookup_model_field_name(field_name), annotation, declared_internal_fields[existing_declared_field_index].second.size(), decomp_info->get_decomp_name(), decomp_info->get_num_chunks());
	}
	
	comp_field_insts.clear();
	for (int i = 0; i < external_procedures_inst->API_specified_field_insts.size(); i ++)
		if (words_are_the_same(lookup_model_field_name(field_name), external_procedures_inst->API_specified_field_insts[i]->get_field_name()))
			comp_field_insts.push_back(external_procedures_inst->API_specified_field_insts[i]);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, comp_field_insts.size() <= 1, "Software error in External_procedure_inst::declare_para_field");
	if ((type_inout == CCPL_PARA_TYPE_IN || type_inout == CCPL_PARA_TYPE_INOUT) && comp_field_insts.size() == 0) 
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the external procedure wants to import this field from the component model, while no instance of this field has been specified when initializing the corresponding external procedures instance \"%s\" at the model code with the annotation \"%s\". Please vefify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, external_procedures_inst->instance_name, annotation_mgr->get_annotation(external_procedures_inst->instance_id,"registering an instance of external procedures"));
	if (type_inout == CCPL_PARA_TYPE_OUT && comp_field_insts.size() == 0)
		EXECUTION_REPORT(REPORT_WARNING, host_comp_id, false, "WARING happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the external procedure can export this field to the component model, while no instance of this field has been specified when initializing the corresponding external procedures instance \"%s\" at the model code with the annotation \"%s\" (the changes on this field will not be used by the component model). Please confirm.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, external_procedures_inst->instance_name, annotation_mgr->get_annotation(external_procedures_inst->instance_id,"registering an instance of external procedures"));

	comp_or_grid_id = grid_id == -1? host_comp_id : grid_id;
	if (existing_declared_field_index == -1)
		num_existing_chunks = 0;
	else num_existing_chunks = declared_internal_fields[existing_declared_field_index].second.size();
	if (data_size == -1 && comp_field_insts.size() > 0 && comp_field_insts[0]->get_decomp_id() == decomp_id && comp_field_insts[0]->get_grid_id() == grid_id && words_are_the_same(comp_field_insts[0]->get_data_type(), data_type) && are_field_units_the_same(comp_field_insts[0]->get_unit(), field_unit)) {
		use_the_same_field_inst = true;
		declared_internal_field.first = comp_field_insts[0];
	}	
	else {
		use_the_same_field_inst = false;
		declared_internal_field.first = memory_manager->alloc_mem(lookup_model_field_name(field_name), decomp_id, comp_or_grid_id, -instance_id, data_type, field_unit, annotation, false);
		if (num_existing_chunks == 0 && data_size > -1)
			declared_internal_field.first->change_to_registered_without_data_buffers();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (declared_internal_field.first->get_is_registered_model_buf() && data_size != -1) || (!declared_internal_field.first->get_is_registered_model_buf() && data_size == -1), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": all data buffers for this field should be allocated by the model (or C-Coupler) at the same time. Please verify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation);					
	}
	if (comp_field_insts.size() > 0 && (comp_field_insts[0]->get_decomp_id() != decomp_id || comp_field_insts[0]->get_grid_id() != grid_id))
		require_coupling_parameters = true;
	if (decomp_info != NULL && decomp_info->get_num_chunks() > 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, specified_size == declared_internal_field.first->get_chunk_data_buf_size(num_existing_chunks), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the data size (%d) specified by \"num_dims\" and \"dims_size\" is different from the right data size (%d) of the number %d chunk. Please verify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, specified_size, declared_internal_field.first->get_chunk_data_buf_size(num_existing_chunks), num_existing_chunks+1);
	}	
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, specified_size == declared_internal_field.first->get_size_of_field(), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the data size (%d) specified by \"num_dims\" and \"dims_size\" is different from the right data size (%d) of the field. Please verify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, specified_size, declared_internal_field.first->get_size_of_field());
	if (data_size == -1) {
		if (decomp_info != NULL && decomp_info->get_num_chunks() > 0) {
			*data_ptr = declared_internal_field.first->get_chunk_buf(num_existing_chunks);			
			EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "When calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": use number %d chunk (%lx) from the field %x of C-Coupler", procedure_name, field_name, lookup_model_field_name(field_name), annotation, num_existing_chunks, *data_ptr, declared_internal_field.first->get_field_instance_id());
		}
		else *data_ptr = declared_internal_field.first->get_data_buf();
	}
	else {		
		if (decomp_info != NULL && decomp_info->get_num_chunks() > 0) {
			declared_internal_field.first->set_data_buf_from_model(*data_ptr, num_existing_chunks);
			EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "When calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": set model buffer to number %d chunk of the field %x in C-Coupler", procedure_name, field_name, lookup_model_field_name(field_name), annotation, num_existing_chunks, declared_internal_field.first->get_field_instance_id());
		}
		else declared_internal_field.first->set_data_buf_from_model(*data_ptr, -1);
	}

	if (comp_field_insts.size() > 0) {
		if (comp_field_insts[0]->get_grid_id() == -1 && declared_internal_field.first->get_grid_id() != -1)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the corresponding field instance has been registered as a scalar variable (not on a grid), while it is now declared as a variable on the grid \"%s\" in this API call. Please vefify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, declared_internal_field.first->get_grid_name());
		if (comp_field_insts[0]->get_grid_id() != -1 && declared_internal_field.first->get_grid_id() == -1)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the corresponding field instance has been registered as a variable on the grid \"%s\", while it is now declared as scalar variable (not on a grid) in this API call. Please vefify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, comp_field_insts[0]->get_grid_name());
		if (comp_field_insts[0]->get_grid_id() != -1)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid_mgr->search_grid_info(comp_field_insts[0]->get_grid_id())->get_original_CoR_grid()->have_the_same_dimensions_with(original_grid_mgr->search_grid_info(declared_internal_field.first->get_grid_id())->get_original_CoR_grid()), "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the grid (\"%s\") of the corresponding field instance when calling the API \"CCPL_external_procedures_inst_init\" to initialize the instance (\"%s\") of the external procedures (\"%s\") and the grid (\"%s\") for declaring this field have different dimensions. Please vefify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, comp_field_insts[0]->get_grid_name(), external_procedures_inst->instance_name, external_procedures_inst->procedures_name, declared_internal_field.first->get_grid_name());
	}

	if (!use_the_same_field_inst && existing_declared_field_index == -1) {
		if (type_inout == CCPL_PARA_TYPE_IN || type_inout == CCPL_PARA_TYPE_INOUT) {
			model_export_field_insts.push_back(comp_field_insts[0]);
			procedure_import_field_insts.push_back(declared_internal_field.first);
		}
		if ((type_inout == CCPL_PARA_TYPE_OUT || type_inout == CCPL_PARA_TYPE_INOUT) && comp_field_insts.size() > 0){
			model_import_field_insts.push_back(comp_field_insts[0]);
			procedure_export_field_insts.push_back(declared_internal_field.first);
		}
	}

	if (existing_declared_field_index == -1) {
		existing_declared_field_index = declared_internal_fields.size();
		declared_internal_fields.push_back(declared_internal_field);
	}
	for (int i = 0; i < given_num_dims; i ++)
		vec_dims_size.push_back(dims_size[i]);
	if (declared_internal_fields[existing_declared_field_index].second.size() > 0)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, declared_internal_fields[existing_declared_field_index].second[0].size() == given_num_dims, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping) at the model code with the annotation \"%s\": the pointers of the number 1 and %d chunks have different dimension numbers (%d vs %d) while the dimension number must be the same. Please verify.", procedure_name, field_name, lookup_model_field_name(field_name), annotation, declared_internal_fields[existing_declared_field_index].second.size()+1, declared_internal_fields[existing_declared_field_index].second[0].size(), given_num_dims);
	declared_internal_fields[existing_declared_field_index].second.push_back(vec_dims_size);

	if (declared_internal_field.second.size() > 1)
		for (int i = 0; i < declared_internal_field.second.size()-1; i ++)
			if (declared_internal_field.first->get_chunk_buf(i) != NULL || *data_ptr != NULL)
				EXECUTION_REPORT(REPORT_ERROR, host_comp_id, declared_internal_field.first->get_chunk_buf(i) != *data_ptr, "When calling the API \"CCPL_external_procedures_para_declare_field\" in the external procedure \"%s\" for the field \"%s\" (model field \"%s\" after name mapping): the No. %d and No. %d chunks correspond to the same data buffer. Please verify the model code with the annotation \"%s\".", procedure_name, field_name, lookup_model_field_name(field_name), i+1, declared_internal_field.second.size(), annotation);

	return declared_internal_field.first->get_field_instance_id();
}


void External_procedure_inst::get_field_pointer(const char *field_name, void **data_ptr, int is_associated, int num_dims, int *dims_size, int chunk_index, const char *data_type, const char *annotation)
{
	int existing_field_index = -1;
	for (int i = 0; i < declared_internal_fields.size(); i ++)
		if (words_are_the_same(declared_internal_fields[i].first->get_field_name(), lookup_model_field_name(field_name))) {
			existing_field_index = i;
			break;
		}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, existing_field_index != -1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_pointer\" to get the data buffer pointer of the field \"%s\" (model field \"%s\" after name mapping) from the external procedure \"%s\" in the dynamic-link library \"%s\" at the model code with the annotation \"%s\": this external procedure has not declared \"%s\" (model field \"%s\" after name mapping) as a parameter field. Please verify.", field_name, lookup_model_field_name(field_name), procedure_name, dl_name, annotation, field_name, lookup_model_field_name(field_name));
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, chunk_index >= 1 && chunk_index <= declared_internal_fields[existing_field_index].second.size(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_pointer\" to get the data buffer pointer of the field \"%s\" (model field \"%s\" after name mapping) from the external procedure \"%s\" in the dynamic-link library \"%s\" at the model code with the annotation \"%s\": the input parameter \"chunk_index\" (%d) is not between 1 and %d. Please verify.", field_name, lookup_model_field_name(field_name), procedure_name, dl_name, annotation, chunk_index, declared_internal_fields[existing_field_index].second.size());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, declared_internal_fields[existing_field_index].second[0].size() == num_dims, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_pointer\" to get the data buffer pointer of the field \"%s\" (model field \"%s\" after name mapping) from the external procedure \"%s\" in the dynamic-link library \"%s\" at the model code with the annotation \"%s\": the parameter \"data_pointer\" of this API call is %d-D while the data buffers of this field have been declared as %d-D before. Please verify.", field_name, lookup_model_field_name(field_name), procedure_name, dl_name, annotation, num_dims, declared_internal_fields[existing_field_index].second[0].size());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, is_associated == 0, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_pointer\" to get the data buffer pointer of the field \"%s\" (model field \"%s\" after name mapping) from the external procedure \"%s\" in the dynamic-link library \"%s\" at the model code with the annotation \"%s\": the parameter \"data_pointer\" has been associated to a data buffer while it should be an empty pointer. Please verify.", field_name, lookup_model_field_name(field_name), procedure_name, dl_name, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type, declared_internal_fields[existing_field_index].first->get_data_type()), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_pointer\" to get the data buffer pointer of the field \"%s\" (model field \"%s\" after name mapping) from the external procedure \"%s\" in the dynamic-link library \"%s\" at the model code with the annotation \"%s\": the implicit data type (%s) of the given \"data_pointer\" of the API call is different from the implicit data type (%s) when declaring this field. Please verify.", field_name, lookup_model_field_name(field_name), procedure_name, dl_name, annotation, data_type, declared_internal_fields[existing_field_index].first->get_data_type());
	if (declared_internal_fields[existing_field_index].first->get_num_chunks() == 0)
		*data_ptr = declared_internal_fields[existing_field_index].first->get_data_buf();
	else *data_ptr = declared_internal_fields[existing_field_index].first->get_chunk_buf(chunk_index-1);
	for (int i = 0; i < num_dims; i ++)
		dims_size[i] = (declared_internal_fields[existing_field_index].second[chunk_index-1])[i];
}


const char *External_procedure_inst::lookup_model_field_name(const char *procedure_field_name)
{
	if (fields_names_mapper == NULL)
		return procedure_field_name;
	
	for (int i = 0; i < (*fields_names_mapper).size(); i ++)
		if (words_are_the_same((*fields_names_mapper)[i].second, procedure_field_name))
			return (*fields_names_mapper)[i].first;

	return procedure_field_name;
}



void External_procedure_inst::run(int chunk_index, const char *annotation)
{
	if (inout_interface_before_run != NULL)
		inout_interface_before_run->execute(true, API_ID_EXTERNAL_PROC_INST_RUN, field_update_status, procedure_import_field_insts.size()+procedure_export_field_insts.size()+1, annotation);

	procedure_run(&instance_id,&chunk_index);

	if (inout_interface_after_run != NULL)
		inout_interface_after_run->execute(true, API_ID_EXTERNAL_PROC_INST_RUN, field_update_status, procedure_import_field_insts.size()+procedure_export_field_insts.size()+1, annotation);
}


void External_procedure_inst::get_procedures_import_field_insts(std::vector<Field_mem_info*> &procedures_import_field_insts)
{
	bool has_the_same_field;


	for (int j = 0; j < procedure_import_field_insts.size(); j ++) {
		has_the_same_field = false;
		for (int i = 0; i < procedures_import_field_insts.size(); i ++)
			if (words_are_the_same(procedures_import_field_insts[i]->get_field_name(), procedure_import_field_insts[j]->get_field_name())) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, procedures_import_field_insts[i] == procedure_import_field_insts[j], "Software error in External_procedure_inst::get_procedures_import_field_insts");
				has_the_same_field = true;
			}
		if (has_the_same_field)
			continue;
		procedures_import_field_insts.push_back(procedure_import_field_insts[j]);
	}
}


void External_procedure_inst::get_procedures_export_field_insts(std::vector<Field_mem_info*> &procedures_export_field_insts)
{
	bool has_the_same_field;


	for (int j = 0; j < procedure_export_field_insts.size(); j ++) {
		has_the_same_field = false;
		for (int i = 0; i < procedures_export_field_insts.size(); i ++)
			if (words_are_the_same(procedures_export_field_insts[i]->get_field_name(), procedure_export_field_insts[j]->get_field_name())) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, procedures_export_field_insts[i] == procedure_export_field_insts[j], "Software error in External_procedure_inst::get_procedures_export_field_insts");
				has_the_same_field = true;
			}
		if (has_the_same_field)
			continue;
		procedures_export_field_insts.push_back(procedure_export_field_insts[j]);
	}
}


void External_procedures_operation::initialize_data(External_procedures_inst *external_procedures_inst, int comp_id, int num_omp_levels)
{
	this->external_procedures_inst = external_procedures_inst;
	this->external_procedure_inst = NULL;
	this->loop_count = 1;
	this->num_total_threads = 0;
	this->num_session_threads = 0;
	this->num_omp_levels = num_omp_levels;
	this->host_comp_id = comp_id;
}


External_procedures_operation::External_procedures_operation(External_procedures_inst *external_procedures_inst, int comp_id, const char *procedure_name, const char *annotation)
{
	initialize_data(external_procedures_inst, comp_id, 0);
	external_procedure_inst = external_procedures_inst->initialize_procedure_inst(procedure_name, external_procedures_inst->dl_name, NULL, NULL, annotation);
}


External_procedures_operation::External_procedures_operation(External_procedures_inst *external_procedures_inst, int comp_id, TiXmlElement *XML_element, const char *XML_file_name, int num_omp_levels, const char *package_fields_name_mapper_name, int package_level, const char *annotation)
{
	initialize_data(external_procedures_inst, comp_id, num_omp_levels);
	
	if (words_are_the_same(XML_element->Value(), "loop") || words_are_the_same(XML_element->Value(), "external_procedures_package")) {
		if (words_are_the_same(XML_element->Value(), "loop")) {
			const char *count_str = get_XML_attribute(comp_id, -1, XML_element, "count", XML_file_name, line_number, "the count of the corresponding loop", "external procedures configuration", true);
			EXECUTION_REPORT(REPORT_ERROR, comp_id, sscanf(count_str, "%d", &loop_count) == 1, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the value (\"%s\") of the attribute \"count\" (loop count) is not an integer. Please verify.", XML_element->Row(), XML_file_name, count_str);
			EXECUTION_REPORT(REPORT_ERROR, comp_id, loop_count >= 1, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the value (\"%d\") of the attribute \"count\" (loop count) should not be smaller than 1. Please verify.", XML_element->Row(), XML_file_name, loop_count);
		}
		else {
			package_fields_name_mapper_name = get_XML_attribute(comp_id, -1, XML_element, "fields_name_mapper", XML_file_name, line_number, "the dynamic-link library of the external procedure", "external procedures configuration", false);
			package_level ++;
		}	
		EXECUTION_REPORT(REPORT_ERROR, comp_id, package_level <= 1, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: an external procedure package is nested in another external procedure package, which is not allowed. Please verify.", XML_element->Row(), XML_file_name);
		for (TiXmlNode *child_node = XML_element->FirstChild(); child_node != NULL; child_node = child_node->NextSibling()) {
			TiXmlElement *child_element = child_node->ToElement();
			if (is_XML_setting_on(comp_id, child_element, XML_file_name, "status of an entry for setting external procedures", "external procedures configuration"))
				children_operations.push_back(new External_procedures_operation(external_procedures_inst, comp_id, child_element, XML_file_name, num_omp_levels, package_fields_name_mapper_name, package_level, annotation));
		}
	}
	else if (words_are_the_same(XML_element->Value(), "openmp_parallel")) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_omp_levels == 0, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the corresponding XML node of \"openmp_parallel\" is nested in another XML node of \"openmp_parallel\" while this kind of nesting is not allowed. Please verify.", XML_element->Row(), XML_file_name);		
		this->num_omp_levels ++;
		for (TiXmlNode *child_node = XML_element->FirstChild(); child_node != NULL; child_node = child_node->NextSibling()) {
			TiXmlElement *child_element = child_node->ToElement();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(child_element->Value(), "parallel_session"), "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the name (\"%s\") of the corresponding XML node is illegal (a right name should be \"parallel_session\" because the XML node of \"openmp_parallel\" can only have the children of \"parallel_session\". Please verify.", child_element->Row(), XML_file_name, child_element->Value());
			if (is_XML_setting_on(comp_id, child_element, XML_file_name, "status of a \"parallel session\"", "external procedures configuration"))
				children_operations.push_back(new External_procedures_operation(external_procedures_inst, comp_id, child_element, XML_file_name, this->num_omp_levels, package_fields_name_mapper_name, package_level, annotation));
		}		
		for (int i = 0; i < children_operations.size(); i ++)
			this->num_total_threads += children_operations[i]->get_num_session_threads();
	}
	else if (words_are_the_same(XML_element->Value(), "parallel_session")) {
		const char *num_threads_str = get_XML_attribute(comp_id, -1, XML_element, "num_threads", XML_file_name, line_number, "the number of threads for the corresponding OpenMP parallel session", "external procedures configuration", true);
		EXECUTION_REPORT(REPORT_ERROR, comp_id, sscanf(num_threads_str, "%d", &this->num_session_threads) == 1, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the value (\"%s\") of the attribute \"num_threads\" (the number of threads for the corresponding OpenMP parallel session) is not an integer. Please verify.", XML_element->Row(), XML_file_name, num_threads_str);
		EXECUTION_REPORT(REPORT_ERROR, comp_id, num_session_threads >= 1, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the value (\"%d\") of the attribute \"num_threads\" (the number of threads for the corresponding OpenMP parallel session) should not be smaller than 1. Please verify.", XML_element->Row(), XML_file_name, loop_count);
		for (TiXmlNode *child_node = XML_element->FirstChild(); child_node != NULL; child_node = child_node->NextSibling()) {
			TiXmlElement *child_element = child_node->ToElement();
			if (is_XML_setting_on(comp_id, child_element, XML_file_name, "status of an entry for setting external procedures", "external procedures configuration"))
				children_operations.push_back(new External_procedures_operation(external_procedures_inst, comp_id, child_element, XML_file_name, num_omp_levels, package_fields_name_mapper_name, package_level, annotation));
		}		
	}
	else if (words_are_the_same(XML_element->Value(), "procedure")) {
		const char *procedure_name = get_XML_attribute(comp_id, -1, XML_element, "name", XML_file_name, line_number, "the name of the external procedure", "external procedures configuration", true);
		const char *dl_name = get_XML_attribute(comp_id, -1, XML_element, "library", XML_file_name, line_number, "the dynamic-link library of the external procedure", "external procedures configuration", false);
		const char *fields_name_mapper_name = get_XML_attribute(comp_id, -1, XML_element, "fields_name_mapper", XML_file_name, line_number, "the dynamic-link library of the external procedure", "external procedures configuration", false);
		if (fields_name_mapper_name == NULL)
			fields_name_mapper_name = package_fields_name_mapper_name;
		std::vector<std::pair<const char *, const char*> > * fields_name_mapper = NULL;
		if (dl_name == NULL)
			dl_name = external_procedures_inst->dl_name;
		if (fields_name_mapper_name != NULL) {
			fields_name_mapper = external_procedures_mgr->search_fields_name_mapper(XML_file_name, fields_name_mapper_name);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, fields_name_mapper != NULL, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures package \"%s\" at the model code with the annotation \"%s\" based on the XML file \"%s\": the fields name mapper \"%s\" for the external procedure \"%s\" is unknown. Please verify.", external_procedures_inst->instance_name, external_procedures_inst->procedures_name, annotation, XML_file_name, fields_name_mapper_name, procedure_name);
		}
		external_procedure_inst = external_procedures_inst->initialize_procedure_inst(procedure_name, dl_name, XML_file_name, fields_name_mapper, annotation);
	}
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Detect an ERROR arround the line number %d in the XML configuration file \"%s\" for external procedures configuration: the name (\"%s\") of the corresponding XML node is illegal (a right name should be \"loop\", \"openmp_parallel\", \"parallel_session\", or \"procedure\"). Please verify.", XML_element->Row(), XML_file_name, XML_element->Value());
}


External_procedures_operation::~External_procedures_operation()
{
	for (int i = 0; i < children_operations.size(); i ++)
		delete children_operations[i];
}


void External_procedures_operation::run(int chunk_index, const char * annotation)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, loop_count >= 1 && !(children_operations.size() > 0 && external_procedure_inst != NULL), "Software error in External_procedures_operation::run");
	for (int j = 0; j < loop_count; j ++) {
		for (int i = 0; i < children_operations.size(); i ++)
			children_operations[i]->run(chunk_index, annotation);
		if (external_procedure_inst != NULL)
			external_procedure_inst->run(chunk_index, annotation);
	}
}


External_procedures_inst::External_procedures_inst(int instance_id, const char *inst_name, const char *procedures_name, const char *ptype, int host_comp_id, const char *dl_name, int process_active, int size_controls, int size_decomps, int size_grids, int size_field_inst, int size_timers, 
                                                   const int *control_vars, const int *decomp_ids, const int *grid_ids, const int *field_inst_ids, const int *timer_ids, const char *annotation)
{
	int tmp_int, num_total_active_procs;


	this->finalized = false;
	check_and_verify_name_format_of_string_for_API(host_comp_id, inst_name, API_ID_EXTERNAL_PROC_INST_INIT, "the name of the external procedure(s) instance", annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(ptype,"individual") || words_are_the_same(ptype,"package"), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the given value (\"%s\") of the input parameter \"%s\" is not \"individual\" and \"package\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, ptype, annotation);	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, strlen(procedures_name) > 0, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures: the name of the external procedures should not be empty. Please verify at the model code with the annotation \"%s\".", inst_name, annotation);
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "External_procedures_inst::External_procedures_inst");
    synchronize_comp_processes_for_API(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", annotation);
    check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", inst_name, "inst_name", annotation);
    check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", procedures_name, "procedures_name", annotation);
    check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", ptype, "ptype", annotation);
	tmp_int = strlen(dl_name) == 0? 0: 1;
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", tmp_int, "dl_name", annotation);	
	if (tmp_int > 0)
		check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", dl_name, "dl_name", annotation);
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_controls, "array size of \"control_vars\" (array size)", annotation);	
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_field_inst, "array size of \"field_inst_ids\" (array size)", annotation);
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_decomps, "array size of \"decomp_ids\" (array size)", annotation);
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_grids, "array size of \"grid_ids\" (array size)", annotation);
	check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_timers, "array size of \"timer_ids\" (array size)", annotation);	
	check_API_parameter_data_array(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", size_controls, sizeof(int), (const char*) control_vars, "control_vars", annotation);
	for (int i = 0; i < size_controls; i ++)
		this->control_vars.push_back(control_vars[i]);
    for (int i = 0; i < size_field_inst; i ++) {
        check_API_parameter_field_instance(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", field_inst_ids[i], "field_inst_ids", annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, memory_manager->get_field_instance(field_inst_ids[i])->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d field instance specified by the input parameter \"field_ids\"\" does not belong to the component model corresponding to the input parameter of \"host_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
		this->API_specified_field_insts.push_back(memory_manager->get_field_instance(field_inst_ids[i]));
    }
	for (int i = 0; i < size_timers; i ++) {
        check_API_parameter_timer(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", timer_ids[i], "timer_ids", annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, timer_mgr->get_timer(timer_ids[i])->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d timer specified by the input parameter \"timer_ids\"\" does not belong to the component model corresponding to the input parameter of \"host_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
		this->timer_ids.push_back(timer_ids[i]);
	}	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, size_grids == size_decomps, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the array size of the input parameter \"grid_ids\" (%d) is different from the array size of the input parameter \"decomp_ids\" (%d). Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, size_grids, size_decomps, annotation);
	for (int i = 0; i < size_grids; i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid_mgr->is_grid_id_legal(grid_ids[i]), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d element of the input parameter \"grid_ids\" is not a legal ID of grids Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid_mgr->search_grid_info(grid_ids[i])->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d grid specified by the input parameter \"grid_ids\"\" does not belong to the component model corresponding to the input parameter of \"host_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);				
		check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", original_grid_mgr->search_grid_info(grid_ids[i])->get_grid_name(), "\"grid_ids\"", annotation);
		if (decomp_ids[i] == -1 || decomp_ids[i] == CCPL_NULL_INT)
			tmp_int = -1;
		else tmp_int = 0;
		check_API_parameter_int(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", tmp_int, "\"decomp_ids\"", annotation);	
		if (decomp_ids[i] != -1 && decomp_ids[i] != CCPL_NULL_INT) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, decomps_info_mgr->is_decomp_id_legal(decomp_ids[i]), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d element of the input parameter \"decomp_ids\" is not a legal ID of parallel decomposition. Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_comp_id() == host_comp_id, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d parallel decomposition specified by the input parameter \"decomp_ids\"\" does not belong to the component model corresponding to the input parameter of \"host_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, original_grid_mgr->search_grid_info(decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_grid_id())->get_original_CoR_grid()->is_subset_of_grid(original_grid_mgr->search_grid_info(grid_ids[i])->get_original_CoR_grid()), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the number %d parallel decomposition in the input parameter \"decomp_ids\"\" is not consistent with the corresponding grid speicified in the input parameter \"grid_ids\". Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, i+1, annotation);
			check_API_parameter_string(host_comp_id, API_ID_EXTERNAL_PROC_INST_INIT, comm, "registering an instance of external procedures", decomps_info_mgr->get_decomp_info(decomp_ids[i])->get_decomp_name(), "\"decomp_ids\"", annotation);
			
		}
		this->grid_ids.push_back(grid_ids[i]);
		this->decomp_ids.push_back(decomp_ids[i]);
	}

	this->host_comp_id = host_comp_id;
	this->instance_name = strdup(inst_name);
	this->procedures_name = strdup(procedures_name);
	this->instance_id = instance_id;
	this->current_process_active = (process_active == 1);
	if (words_are_the_same(ptype,"individual"))
		this->procedures_type = 0;
	else this->procedures_type = 1;	
	MPI_Reduce(&process_active, &num_total_active_procs, 1, MPI_INT, MPI_SUM, 0, comm);
	if (comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_current_proc_local_id() == 0)
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, num_total_active_procs > 0, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": no process is active to run the instance of external procedures. Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, annotation);
	EXECUTION_REPORT(REPORT_ERROR,-1, MPI_Comm_split(comm, process_active, 0, &(this->local_comm)) == MPI_SUCCESS);
	if (!this->current_process_active) {
		this->local_comm = MPI_COMM_NULL;
		this->num_proc_in_local_comm = 0;
		this->proc_id_in_local_comm = -1;
	}
	else {
		MPI_Comm_rank(this->local_comm, &proc_id_in_local_comm);
		MPI_Comm_size(this->local_comm, &num_proc_in_local_comm);
	}

	annotation_mgr->add_annotation(instance_id, "registering an instance of external procedures", annotation);
	external_procedures_mgr->add_external_procedures_inst(this);

	for (int i = 0; i < API_specified_field_insts.size(); i ++)
		for (int j = i+1; j < API_specified_field_insts.size(); j ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !words_are_the_same(API_specified_field_insts[i]->get_field_name(),API_specified_field_insts[j]->get_field_name()), "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the input parameter \"field_inst_ids\" contains multiple elements corresponding to the same field \"%s\" (at number %d and %d elements). Please verify at the model code with the annotation \"%s\" (please specify at most one element for this field).", inst_name, procedures_name, API_specified_field_insts[i]->get_field_name(), i+1, j+1, annotation);

	this->dl_name = NULL;
	this->XML_file_name = NULL;
	if (this->procedures_type == 0) {
		if (strlen(dl_name) != 0)
			this->dl_name = strdup(dl_name);
		root_operation = new External_procedures_operation(this, host_comp_id, procedures_name, annotation);
	}
	else {
		TiXmlElement *config_XML_element = get_XML_file_with_configuration();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, config_XML_element != NULL, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures package \"%s\" at the model code with the annotation \"%s\": no corresponding configuration information can be found from the corresponding XML files under the directory \"%s\". Please verify.", inst_name, procedures_name, annotation, comp_comm_group_mgt_mgr->get_external_procedure_config_dir());
		const char *config_dl_name = get_XML_attribute(host_comp_id, -1, config_XML_element, "library", XML_file_name, line_number, "the dynamic-link library of the external procedure", "external procedures configuration", false);
		if (config_dl_name != NULL)
			this->dl_name = strdup(config_dl_name);
		else if (strlen(dl_name) != 0)
			this->dl_name = strdup(dl_name);
		root_operation = new External_procedures_operation(this, host_comp_id, config_XML_element, XML_file_name, 0, NULL, 0, annotation);
	}
}


External_procedures_inst::~External_procedures_inst()
{
	delete [] instance_name;
	delete [] procedures_name;
	if (dl_name != NULL)
		delete [] dl_name;
	if (XML_file_name != NULL)
		delete [] XML_file_name;
	if (root_operation != NULL)
		delete root_operation;
}


TiXmlElement *External_procedures_inst::get_XML_file_with_configuration()
{
	char current_XML_file_name[NAME_STR_SIZE];
	Comp_comm_group_mgt_node *host_comp_node = comp_comm_group_mgt_mgr->search_global_node(host_comp_id);
	TiXmlDocument *current_XML_file;

	for (Comp_comm_group_mgt_node *current_comp_node = host_comp_node; current_comp_node != NULL; current_comp_node = current_comp_node->get_parent()) {
		if (words_are_the_same(current_comp_node->get_comp_type(), COMP_TYPE_PSEUDO_COUPLED))
			continue;
		if (words_are_the_same(current_comp_node->get_comp_type(), COMP_TYPE_ROOT))
			sprintf(current_XML_file_name, "%s/overall.external_procedures.xml", comp_comm_group_mgt_mgr->get_external_procedure_config_dir());
		else sprintf(current_XML_file_name, "%s/%s.external_procedures.xml", comp_comm_group_mgt_mgr->get_external_procedure_config_dir(), current_comp_node->get_comp_full_name());
		current_XML_file = open_XML_file_to_read(host_comp_id, current_XML_file_name, host_comp_node->get_comm_group(), false);
		if (current_XML_file == NULL)
			continue;
		TiXmlElement *XML_element = current_XML_file->FirstChildElement();
		for (TiXmlNode *XML_node = get_XML_first_child_of_unique_root(host_comp_id, XML_file_name, current_XML_file); XML_node != NULL; XML_node = XML_node->NextSibling()) {
			if (XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
				continue;
			XML_element = XML_node->ToElement();
			if (!words_are_the_same(XML_element->Value(),"external_procedures_package"))
				continue;
			if (!is_XML_setting_on(host_comp_id, XML_element, XML_file_name, "the status \"external_procedures_package\"", "external procedures configuration"))
				continue;			
			const char *current_procedures_name = get_XML_attribute(host_comp_id, -1, XML_element, "name", XML_file_name, line_number, "the name of the external procedures package", "external procedures configuration", true);
			if (words_are_the_same(current_procedures_name, procedures_name)) {
				XML_file_name = strdup(current_XML_file_name);
				XML_file = current_XML_file;
				EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "In the process of calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures package \"%s\": find the XML file \"%s\" with the configuration information of the external procedures package", instance_name, procedures_name, XML_file_name);
				return XML_element;
			}
		}
		delete current_XML_file;
	}

	return NULL;
}


int External_procedures_inst::declare_para_field(int procedure_inst_id, void **data_ptr, const char *field_name, const char *data_type, int type_inout, int decomp_id, int grid_id, int *dims_size, int size_dim_array, int pointer_num_dims, int data_size, const char *field_unit, const char *annotation)
{	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, type_inout == CCPL_PARA_TYPE_IN || type_inout == CCPL_PARA_TYPE_OUT || type_inout == CCPL_PARA_TYPE_INOUT, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\": the value (currently is %d) of the input parameter \"type_inout\" must be CCPL_PARA_TYPE_IN (1), CCPL_PARA_TYPE_OUT (2) or CCPL_PARA_TYPE_INOUT (3). Please verify the model code with the annotation of \"%s\"" , field_name, type_inout, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, size_dim_array >= pointer_num_dims, "ERROR happens when calling the API \"CCPL_external_procedures_para_declare_field\" for the field \"%s\": the array size (%d) of the input parameter \"dims_size\" is smaller than the number of dimensions of the given pointer (%d). Please verify the model code with the annotation of \"%s\"" , field_name, size_dim_array, pointer_num_dims, annotation);
	
	return external_procedure_insts[GET_INST_PROCEDURE_INDEX(procedure_inst_id)]->declare_para_field(data_ptr, field_name, data_type, type_inout, decomp_id, grid_id, pointer_num_dims, dims_size, data_size, field_unit, annotation);
}


External_procedure_inst *External_procedures_inst::initialize_procedure_inst(const char *procedure_name, const char *dl_name, const char *XML_file_name, std::vector<std::pair<const char *, const char*> > *fields_name_mapper, const char *annotation)
{
	if (XML_file_name != NULL)
		EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Start to initialize a instance of the external procedure \"%s\" based on an XML configuration file \"%s\"", procedure_name, XML_file_name);
	return new External_procedure_inst(this, instance_id | (external_procedure_insts.size() << 12), instance_name, procedure_name, host_comp_id, dl_name, fields_name_mapper, annotation);
}


void External_procedures_inst::get_field_pointer(int procedure_inst_id, const char *field_name, void **data_ptr, int is_associated, int num_dims, int *dims_size, int chunk_index, const char *data_type, const char *annotation)
{
	external_procedure_insts[GET_INST_PROCEDURE_INDEX(procedure_inst_id)]->get_field_pointer(field_name, data_ptr, is_associated, num_dims, dims_size, chunk_index, data_type, annotation);
}


void External_procedures_inst::run(int chunk_index, const char *annotation)
{
	if (!current_process_active)
		return;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, chunk_index == -1 || chunk_index > 0, "ERROR happens when calling the API \"CCPL_external_procedures_inst_run\" to run the instance (\"%s\") of the external procedures \"%s\": the value (%d) of \"chunk_index\" is not -1 or a possitive value. Please verify the model code with the annotation of \"%s\".", instance_name, procedures_name, chunk_index, annotation);
	if (finalized)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !finalized, "ERROR happens when calling the API \"CCPL_external_procedures_inst_run\" to run the instance (\"%s\") of the external procedures \"%s\" at the model code with the annotation \"%s\": this instance of external procedures has been finalized before at the model code with the annotation \"%s\". Please verify", instance_name, procedures_name, annotation, annotation_mgr->get_annotation(instance_id, "finalization"));
	if (chunk_index > 0)
		for (int i = 0; i < external_procedure_insts.size(); i ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !external_procedure_insts[i]->does_require_coupling_parameters(), "ERROR happens when calling the API \"CCPL_external_procedures_inst_run\" to run the instance (\"%s\") of the external procedures \"%s\": the current value (%d) of \"chunk_index\" indicates that only one chunk of the external procedures will be executed, but coupling for the declared fields will be required, which requires to execute all chunks at the same time (the value of \"chunk_index\" should be -1). Please verify the model code with the annotation of \"%s\".", instance_name, procedures_name, chunk_index, annotation);
	root_operation->run(chunk_index, annotation);
}


void External_procedures_inst::finalize(const char *annotation)
{
	if (current_process_active)
	    synchronize_comp_processes_for_API(host_comp_id, API_ID_EXTERNAL_PROC_INST_FINALIZE, local_comm, "finalizing an instance of external procedures", annotation);
	if (finalized)
		return;

	annotation_mgr->add_annotation(instance_id, "finalization", annotation);
	finalized = true;
	delete root_operation;
	root_operation = NULL;
	for (int i = 0; i < external_procedure_insts.size(); i ++)
		delete external_procedure_insts[i];
}


void External_procedures_inst::get_procedures_import_field_insts(std::vector<Field_mem_info*> &procedures_import_field_insts)
{
	procedures_import_field_insts.clear();
	for (int i = 0; i < external_procedure_insts.size(); i ++)
		external_procedure_insts[i]->get_procedures_import_field_insts(procedures_import_field_insts);
}


void External_procedures_inst::get_procedures_export_field_insts(std::vector<Field_mem_info*> &procedures_export_field_insts)
{
	procedures_export_field_insts.clear();
	for (int i = 0; i < external_procedure_insts.size(); i ++)
		external_procedure_insts[i]->get_procedures_export_field_insts(procedures_export_field_insts);
}


External_procedures_mgt::External_procedures_mgt()
{
}


External_procedures_mgt::~External_procedures_mgt()
{
	for (int i = 0; i < dl_handlers_info.size(); i ++) {
		delete [] dl_handlers_info[i].second;
		dlclose(dl_handlers_info[i].first);
	}
	
	for (int i = 0; i < registered_external_procedures_insts.size(); i ++)
		delete registered_external_procedures_insts[i];

	for (int i = 0; i < models_external_procedure_fields_names_mappers.size(); i ++)
		delete models_external_procedure_fields_names_mappers[i];
}


void *External_procedures_mgt::get_or_open_dl_handler(int host_comp_id, const char *inst_name, const char *procedures_name, const char *dl_name, char **final_dl_full_name, const char *annotation)
{
	char tmp_string[NAME_STR_SIZE];
	std::pair<void*, char*> dl_handler_info;

	if (dl_name != NULL && strlen(dl_name) > 0)
		sprintf(tmp_string, "%s/%s", comp_comm_group_mgt_mgr->get_external_procedure_lib_dir(), dl_name);
	else sprintf(tmp_string, "%s/%s", comp_comm_group_mgt_mgr->get_external_procedure_lib_dir(), "CCPL_public_external_procedures.so");
	*final_dl_full_name = strdup(tmp_string);

	for (int i = 0; i < dl_handlers_info.size(); i ++)
		if (words_are_the_same(dl_handlers_info[i].second, tmp_string))
			return dl_handlers_info[i].first;
		
	dl_handler_info.second = strdup(tmp_string);
	dl_handler_info.first = dlopen(tmp_string, RTLD_NOW);
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, dl_handler_info.first != NULL, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\": the specified dynamic-link library (%s) cannot be opened (%s). Please verify at the model code with the annotation \"%s\".", inst_name, procedures_name, tmp_string, dlerror(), annotation);
	dl_handlers_info.push_back(dl_handler_info);

	return dl_handler_info.first;
}


int External_procedures_mgt::initialize_external_procedures_inst(const char *inst_name, const char *procedures_name, const char *ptype, int host_comp_id, const char *dl_name, int process_active, int size_controls, int size_decomps, int size_grids, int size_field_inst, int size_timers, 
                                                                 const int *control_vars, const int *decomp_ids, const int *grid_ids, const int *field_inst_ids, const int *timer_ids, const char *annotation)
{
	for (int i = 0; i < registered_external_procedures_insts.size(); i ++) {
		if (registered_external_procedures_insts[i]->has_been_finalized())
			continue;
		if (registered_external_procedures_insts[i]->get_host_comp_id() == host_comp_id && words_are_the_same(registered_external_procedures_insts[i]->get_instance_name(), inst_name)) 
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "ERROR happens when calling the API \"CCPL_external_procedures_inst_init\" to initialize an instance \"%s\" of external procedures \"%s\" at the model code with the annotation \"%s\": annother instance with the same name has been initialized before (at the model code with the annotation \"%s\"). Please verify.", inst_name, procedures_name, annotation, annotation_mgr->get_annotation(registered_external_procedures_insts[i]->get_instance_id(), "registering an instance of external procedures")); 
	}

	int instance_id = TYPE_EXT_PROCEDURE_PREFIX | (registered_external_procedures_insts.size());
	External_procedures_inst *procedures_inst = new External_procedures_inst(instance_id, inst_name, procedures_name, ptype, host_comp_id, dl_name, process_active, size_controls, size_decomps, size_grids, size_field_inst, size_timers, 
																				 control_vars, decomp_ids, grid_ids, field_inst_ids, timer_ids, annotation);
	
	return instance_id;
}


External_procedures_inst *External_procedures_mgt::get_procedures_inst(int instance_id, int API_id, const char *annotation)
{
	char API_label[NAME_STR_SIZE]; 
	int instance_index;

	
	get_API_hint(-1, API_id, API_label);
	instance_index = GET_PROCEDURES_INST_INDEX(instance_id);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (instance_id&TYPE_ID_PREFIX_MASK) == TYPE_EXT_PROCEDURE_PREFIX && instance_index < registered_external_procedures_insts.size(), "ERROR happens when calling the API \"%s\": the given \"instance_id\" (%x) is not a legal ID of an instance of external procedures. Please verify the model code with the anotation \"%s\"", API_label, instance_id, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, GET_INST_PROCEDURE_INDEX(instance_id) < registered_external_procedures_insts[instance_index]->get_num_procedures(), "ERROR happens when calling the API \"%s\": the given \"instance_id\" (%x) is not a legal ID of an instance of external procedures. Please verify the model code with the anotation \"%s\"", API_label, instance_id, annotation);

	return registered_external_procedures_insts[instance_index];
}


MPI_Comm External_procedures_mgt::get_instance_local_comm(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_COMM, annotation)->get_local_comm();
}


int External_procedures_mgt::get_instance_comp_id(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_COMP, annotation)->get_host_comp_id();
}


int External_procedures_mgt::get_instance_process_active(int instance_id, const char *annotation)
{
	if (get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_PROC_ACTIVE, annotation)->is_current_process_active())
		return 1;
	
	return 0;
}


int External_procedures_mgt::get_instance_num_grid_decomps(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_NUM_GRID_DECOMPS, annotation)->get_num_grid_decomps();
}


int External_procedures_mgt::get_instance_grid_id(int instance_id, int grid_index, const char *annotation)
{
	External_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_GRID_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), grid_index >= 1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of grid_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), grid_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), grid_index <= procedures_inst->get_num_grid_decomps(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_grid_ID\" regarding the external procedures instance \"%s\": the value of grid_index (%d) cannot be larger than %d (the number of given grids when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), grid_index, procedures_inst->get_num_grid_decomps(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_grid_id(grid_index);
}


int External_procedures_mgt::get_instance_decomp_id(int instance_id, int decomp_index, const char *annotation)
{
	External_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_DECOMP_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), decomp_index >= 1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of decomp_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), decomp_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), decomp_index <= procedures_inst->get_num_grid_decomps(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_decomp_ID\" regarding the external procedures instance \"%s\": the value of decomp_index (%d) cannot be larger than %d (the number of given parallel decompositions when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), decomp_index, procedures_inst->get_num_grid_decomps(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_decomp_id(decomp_index);
}


int External_procedures_mgt::get_instance_num_control_vars(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_NUM_CONTROLS, annotation)->get_num_control_vars();
}


int External_procedures_mgt::get_instance_control_var(int instance_id, int control_var_index, const char *annotation)
{
	External_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_CONTROL_VAR, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), control_var_index >= 1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of control_var_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), control_var_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), control_var_index <= procedures_inst->get_num_control_vars(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_control_var\" regarding the external procedures instance \"%s\": the value of control_var_index (%d) cannot be larger than %d (the number of given control variables when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), control_var_index, procedures_inst->get_num_control_vars(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_control_var(control_var_index);
}


int External_procedures_mgt::get_instance_num_timers(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_NUM_TIMERS, annotation)->get_num_timers();
}


int External_procedures_mgt::get_instance_timer_id(int instance_id, int timer_index, const char *annotation)
{
	External_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_TIMER_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), timer_index >= 1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of timer_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), timer_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), timer_index <= procedures_inst->get_num_timers(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_timer_ID\" regarding the external procedures instance \"%s\": the value of timer_index (%d) cannot be larger than %d (the number of given timers when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), timer_index, procedures_inst->get_num_timers(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_timer_id(timer_index);
}


int External_procedures_mgt::get_instance_num_fields(int instance_id, const char *annotation)
{
	return get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_NUM_FIELDS, annotation)->get_num_specified_field_instances();
}


int External_procedures_mgt::get_instance_field_id(int instance_id, int field_index, const char *annotation)
{
	External_procedures_inst *procedures_inst = get_procedures_inst(instance_id, API_ID_EXTERNAL_PROC_INST_GET_FIELD_ID, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), field_index >= 1, "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of field_index (%d) should be larger than 0. Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), field_index, annotation);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, procedures_inst->get_host_comp_id(), field_index <= procedures_inst->get_num_specified_field_instances(), "ERROR happens when calling the API \"CCPL_external_procedures_para_get_field_ID\" regarding the external procedures instance \"%s\": the value of field_index (%d) cannot be larger than %d (the number of given fields when initializing the instance at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\"", procedures_inst->get_instance_name(), field_index, procedures_inst->get_num_specified_field_instances(), annotation_mgr->get_annotation(procedures_inst->get_instance_id(),"registering an instance of external procedures"), annotation);
	return procedures_inst->get_specified_field_instance_id(field_index);
}


void External_procedures_mgt::add_model_external_procedure_fields_names_mappers(int comp_id)
{
	char XML_file_name[NAME_STR_SIZE];
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	bool has_fields_names_mappers = false;
	TiXmlDocument *XML_file;


	if (comp_id != -1) {
	    sprintf(XML_file_name, "%s/%s.external_procedures.xml", comp_comm_group_mgt_mgr->get_external_procedure_config_dir(), comp_node->get_comp_full_name());
		XML_file = open_XML_file_to_read(comp_id, XML_file_name, comp_node->get_comm_group(), false);
	}
	else {
	    sprintf(XML_file_name, "%s/overall.external_procedures.xml", comp_comm_group_mgt_mgr->get_external_procedure_config_dir());
		XML_file = open_XML_file_to_read(comp_id, XML_file_name, MPI_COMM_WORLD, false);
	}
	if (XML_file == NULL)
		return;

	TiXmlElement *XML_element = XML_file->FirstChildElement();
	for (TiXmlNode *XML_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file); XML_node != NULL; XML_node = XML_node->NextSibling()) {
		if (XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
			continue;
		XML_element = XML_node->ToElement();
		if (!words_are_the_same(XML_element->Value(),"fields_name_mappers"))
			continue;
		if (!is_XML_setting_on(comp_id, XML_element, XML_file_name, "the status \"fields_name_mappers\"", "external procedures configuration"))
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !has_fields_names_mappers, "ERROR happens when loading the XML file \"%s\" for external procedures configuration: this XML file contains at least two active sections of \"fields_name_mappers\", while at most one active section is allowed. Please verify.", XML_file_name);
		models_external_procedure_fields_names_mappers.push_back(new External_procedure_fields_names_mappers(comp_id, XML_element, XML_file_name));		
		has_fields_names_mappers = true;
	}
	
	delete XML_file;
}


std::vector<std::pair<const char *, const char*> > *External_procedures_mgt::search_fields_name_mapper(const char *XML_file_name, const char *mapper_name)
{
	for (int i = 0; i < models_external_procedure_fields_names_mappers.size(); i ++)
		if (words_are_the_same(models_external_procedure_fields_names_mappers[i]->get_XML_file_name(), XML_file_name))
			return models_external_procedure_fields_names_mappers[i]->search_fields_name_mapper(mapper_name);

	return NULL;
}

