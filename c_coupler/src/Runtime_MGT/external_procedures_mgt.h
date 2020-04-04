/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef EXTERNAL_PROCEDURES_MGT
#define EXTERNAL_PROCEDURES_MGT


#include <mpi.h>
#include <vector>
#include "memory_mgt.h"
#include "inout_interface_mgt.h"


#define CCPL_PARA_TYPE_IN		1
#define CCPL_PARA_TYPE_OUT		2
#define CCPL_PARA_TYPE_INOUT	3


#define GET_PROCEDURES_INST_INDEX(ID)       (ID & 0x00000FFF)  
#define GET_INST_PROCEDURE_INDEX(ID)        ((ID & 0x00FFF000) >> 12)


class External_procedures_inst;


class External_procedure_fields_names_mappers
{
	private:
		const char *XML_file_name;
		std::vector<std::pair<const char *, std::vector<std::pair<const char *, const char*> > > >  fields_names_mappers;

	public:
		External_procedure_fields_names_mappers(int, TiXmlElement *, const char *);
		const char *get_XML_file_name() { return XML_file_name; }
		std::vector<std::pair<const char *, const char*> > *search_fields_name_mapper(const char*);
		~External_procedure_fields_names_mappers();
};


class External_procedure_inst
{
	private:
		int instance_id;
		int host_comp_id;
		const char *dl_name;
		void *local_dl_handler;
		const char *procedure_name;
		void (*procedure_init)(int *);
		void (*procedure_run)(int *, int*);
		void (*procedure_finalize)(int *);
		External_procedures_inst *external_procedures_inst;
		std::vector<Field_mem_info*> model_import_field_insts;
		std::vector<Field_mem_info*> model_export_field_insts;
		std::vector<Field_mem_info*> procedure_import_field_insts;
		std::vector<Field_mem_info*> procedure_export_field_insts;
		std::vector<std::pair<Field_mem_info*, std::vector<std::vector<int> > > > declared_internal_fields;
		Inout_interface *inout_interface_before_run;
		Inout_interface *inout_interface_after_run;
		int *field_update_status;
		std::vector<std::pair<const char *, const char*> > *fields_names_mapper;
		bool require_coupling_parameters;

	public:		
		External_procedure_inst(External_procedures_inst*, int, const char *, const char *, int, const char *, std::vector<std::pair<const char *, const char*> > *, const char *);
		int declare_para_field(void **, const char *, const char *, int, int, int, int, int *, int, const char *, const char *);
		~External_procedure_inst();
		int get_instance_id() { return instance_id; }
		void get_field_pointer(const char *, void **, int, int, int *, int, const char *, const char *);
		const char *lookup_model_field_name(const char*);
		bool does_require_coupling_parameters() { return require_coupling_parameters; }
		void run(int, const char*);
		void get_procedures_import_field_insts(std::vector<Field_mem_info*> &);
		void get_procedures_export_field_insts(std::vector<Field_mem_info*> &);
};


class External_procedures_operation
{
	private:
		int host_comp_id;
		External_procedures_inst *external_procedures_inst;
		External_procedure_inst *external_procedure_inst;
		int loop_count;
		int num_total_threads;
		int num_session_threads;
		int num_omp_levels;
		std::vector<External_procedures_operation *> children_operations;

		int get_num_session_threads() { return num_session_threads; }		
		void initialize_data(External_procedures_inst *, int, int);
		
	public:
		External_procedures_operation(External_procedures_inst *, int, TiXmlElement *, const char *, int, const char *, int, const char*);
		External_procedures_operation(External_procedures_inst *, int, const char *, const char *);
		~External_procedures_operation();
		void run(int, const char *);
};


class External_procedures_inst
{
	private:
		friend class External_procedure_inst;
		friend class External_procedures_operation;
		int instance_id;
		int host_comp_id;
		int procedures_type;    // 0: individual procedure; 1: procedure package
		const char *instance_name;
		const char *procedures_name;
		const char *dl_name;
		const char *XML_file_name;
		TiXmlDocument *XML_file;
		bool current_process_active;
		MPI_Comm local_comm;
		int num_proc_in_local_comm;
		int proc_id_in_local_comm;
		std::vector<int> control_vars;
		std::vector<int> decomp_ids;
		std::vector<int> grid_ids;
		std::vector<Field_mem_info*> API_specified_field_insts;
		std::vector<int> timer_ids;
		std::vector<External_procedure_inst*> external_procedure_insts;
		External_procedures_operation *root_operation;
		bool finalized;

	public:		
		External_procedures_inst(int, const char *, const char *, const char *, int, const char *, int, int, int, int, int, int, const int *, const int *, const int *, const int *, const int *, const char *);
		int get_host_comp_id() { return host_comp_id; }
		const char *get_instance_name() { return instance_name; }
		const char *get_procedures_name() { return procedures_name; }
		~External_procedures_inst();
		int get_num_grid_decomps() { return decomp_ids.size(); }
		int get_num_control_vars() { return control_vars.size(); }
		int get_num_specified_field_instances() { return API_specified_field_insts.size(); }
		int get_num_timers() { return timer_ids.size(); }
		int get_instance_id() { return instance_id; }
		MPI_Comm get_local_comm() { return local_comm; }
		int get_grid_id(int grid_index) { return grid_ids[grid_index-1]; }
		int get_decomp_id(int decomp_index) { return decomp_ids[decomp_index-1]; }
		int get_control_var(int control_var_index) { return control_vars[control_var_index-1]; }
		int get_timer_id(int timer_index) { return timer_ids[timer_index-1]; }
		int get_specified_field_instance_id(int field_index) { return API_specified_field_insts[field_index-1]->get_field_instance_id(); }
		bool is_current_process_active() { return current_process_active; }
		int declare_para_field(int, void **, const char *, const char *, int, int, int, int *, int, int, int, const char *, const char *);
		void get_field_pointer(int, const char*, void **, int, int, int *, int, const char*, const char*);
		int get_num_procedures() { return external_procedure_insts.size(); }
		void add_external_procedure_inst(External_procedure_inst *instance) { external_procedure_insts.push_back(instance); }
		External_procedure_inst *initialize_procedure_inst(const char *, const char *, const char *, std::vector<std::pair<const char *, const char*> > *, const char *);
		TiXmlElement *get_XML_file_with_configuration();
		void run(int, const char*);
		void finalize(const char*);
		bool has_been_finalized() { return finalized; }
		void get_procedures_import_field_insts(std::vector<Field_mem_info*> &);
		void get_procedures_export_field_insts(std::vector<Field_mem_info*> &);
};


class External_procedures_mgt
{
	private:
		std::vector<External_procedures_inst *> registered_external_procedures_insts;
		std::vector<std::pair<void *, char *> > dl_handlers_info;
		std::vector<External_procedure_fields_names_mappers *> models_external_procedure_fields_names_mappers;

	public:
		External_procedures_mgt();
		~External_procedures_mgt();
		int initialize_external_procedures_inst(const char *, const char *, const char *, int, const char *, int, int, int, int, int, int, const int *, const int *, const int *, const int *, const int *, const char *);
		void *get_or_open_dl_handler(int, const char *, const char *, const char *, char **, const char *);	
		External_procedures_inst *get_procedures_inst(int, int, const char*);
		MPI_Comm get_instance_local_comm(int, const char *);
		int get_instance_comp_id(int, const char *);
		int get_instance_process_active(int, const char *);
		int get_instance_num_grid_decomps(int, const char *);
		int get_instance_grid_id(int, int, const char *);
		int get_instance_decomp_id(int, int, const char *);
		int get_instance_num_control_vars(int, const char *);
		int get_instance_control_var(int, int, const char *);
		int get_instance_num_timers(int, const char *);
		int get_instance_timer_id(int, int, const char *);
		int get_instance_num_fields(int, const char *);
		int get_instance_field_id(int, int, const char *);
		void add_external_procedures_inst(External_procedures_inst *instance) { registered_external_procedures_insts.push_back(instance); }
		void add_model_external_procedure_fields_names_mappers(int);
		std::vector<std::pair<const char *, const char*> > *search_fields_name_mapper(const char *, const char*);
};

#endif

